//! Existing Secret Service over an existing local Unix D-Bus socket only.
//! No autolaunch/unixexec/TCP, Unlock, Prompt, daemon launch or file fallback.
use crate::{Error, MAX_TOKEN_BYTES, Result, SERVICE, Scope, Secret};
use std::{collections::HashMap, time::Duration};
use zbus::{
    blocking::{Connection, Proxy, connection::Builder, proxy::Builder as ProxyBuilder},
    proxy::{CacheProperties, MethodFlags},
    zvariant::{OwnedObjectPath, OwnedValue, Value},
};
const BUS: &str = "org.freedesktop.secrets";
const ROOT: &str = "/org/freedesktop/secrets";
const SERVICE_IFACE: &str = "org.freedesktop.Secret.Service";
type DbusSecret = (OwnedObjectPath, Vec<u8>, Vec<u8>, String);
fn error(error: zbus::Error) -> Error {
    if let zbus::Error::MethodError(name, _, _) = error {
        match name.as_str() {
            "org.freedesktop.Secret.Error.IsLocked" => Error::Locked,
            "org.freedesktop.DBus.Error.AccessDenied" => Error::Denied,
            _ => Error::Unavailable,
        }
    } else {
        Error::Unavailable
    }
}
fn local_address(address: &zbus::Address) -> bool {
    matches!(address.transport(), zbus::address::Transport::Unix(_))
}
struct Bus {
    connection: Connection,
    owner: String,
}
fn unique_owner(owner: &str) -> bool {
    zbus::names::UniqueName::try_from(owner).is_ok()
}
fn connect() -> Result<Bus> {
    let address = zbus::Address::session().map_err(error)?;
    if !local_address(&address) {
        return Err(Error::Unavailable);
    }
    let connection = Builder::address(address)
        .map_err(error)?
        .method_timeout(Duration::from_secs(5))
        .max_queued(16)
        .build()
        .map_err(error)?;
    let bus: Proxy<'_> = ProxyBuilder::new(&connection)
        .destination("org.freedesktop.DBus")
        .map_err(error)?
        .path("/org/freedesktop/DBus")
        .map_err(error)?
        .interface("org.freedesktop.DBus")
        .map_err(error)?
        .cache_properties(CacheProperties::No)
        .build()
        .map_err(error)?;
    // Bind the complete operation to one existing unique owner. It cannot be
    // replaced by another daemon between lookup and a secret-bearing method.
    let owner: String = bus
        .call_with_flags("GetNameOwner", MethodFlags::NoAutoStart.into(), &(BUS,))
        .map_err(error)?
        .ok_or(Error::Unavailable)?;
    if !unique_owner(&owner) {
        return Err(Error::Unavailable);
    }
    drop(bus);
    Ok(Bus { connection, owner })
}
fn proxy<'a>(connection: &Bus, path: &'a str, interface: &'a str) -> Result<Proxy<'a>> {
    ProxyBuilder::new(&connection.connection)
        .destination(connection.owner.clone())
        .map_err(error)?
        .path(path)
        .map_err(error)?
        .interface(interface)
        .map_err(error)?
        // Never let zbus perform an implicit property Get/GetAll. All actual
        // properties use the explicit NoAutoStart call path below.
        .cache_properties(CacheProperties::No)
        .build()
        .map_err(error)
}
fn attrs(scope: &Scope) -> HashMap<&str, &str> {
    HashMap::from([("application", SERVICE), ("account", scope.account())])
}
fn item(connection: &Bus, scope: &Scope) -> Result<OwnedObjectPath> {
    let service = proxy(connection, ROOT, SERVICE_IFACE)?;
    let (unlocked, locked): (Vec<OwnedObjectPath>, Vec<OwnedObjectPath>) = service
        .call_with_flags(
            "SearchItems",
            MethodFlags::NoAutoStart.into(),
            &(attrs(scope),),
        )
        .map_err(error)?
        .ok_or(Error::Unavailable)?;
    if !locked.is_empty() {
        return Err(Error::Locked);
    }
    match unlocked.len() {
        0 => Err(Error::Missing),
        1 => Ok(unlocked.into_iter().next().ok_or(Error::Missing)?),
        _ => Err(Error::Ambiguous),
    }
}
struct Session {
    connection: Bus,
    path: OwnedObjectPath,
}
impl Session {
    fn open(connection: Bus) -> Result<Self> {
        let service = proxy(&connection, ROOT, SERVICE_IFACE)?;
        let (_output, path): (OwnedValue, OwnedObjectPath) = service
            .call_with_flags(
                "OpenSession",
                MethodFlags::NoAutoStart.into(),
                &("plain", Value::from("")),
            )
            .map_err(error)?
            .ok_or(Error::Unavailable)?;
        if path.as_str() == "/" {
            return Err(Error::Unavailable);
        }
        drop(service);
        Ok(Self { connection, path })
    }
}
impl Drop for Session {
    fn drop(&mut self) {
        if let Ok(session) = proxy(
            &self.connection,
            self.path.as_str(),
            "org.freedesktop.Secret.Session",
        ) {
            let _: std::result::Result<Option<()>, _> =
                session.call_with_flags("Close", MethodFlags::NoAutoStart.into(), &());
        }
    }
}
pub fn read(scope: &Scope) -> Result<Secret> {
    let connection = connect()?;
    let path = item(&connection, scope)?;
    let session = Session::open(connection)?;
    let item = proxy(
        &session.connection,
        path.as_str(),
        "org.freedesktop.Secret.Item",
    )?;
    let (returned, parameters, bytes, content): DbusSecret = item
        .call_with_flags(
            "GetSecret",
            MethodFlags::NoAutoStart.into(),
            &(session.path.clone(),),
        )
        .map_err(error)?
        .ok_or(Error::Unavailable)?;
    let bytes = zeroize::Zeroizing::new(bytes);
    if returned != session.path
        || !parameters.is_empty()
        || bytes.len() > MAX_TOKEN_BYTES
        || !matches!(
            content.as_str(),
            "text/plain" | "text/plain; charset=utf-8" | "text/plain;charset=utf-8"
        )
    {
        return Err(Error::InvalidSecret);
    }
    Secret::new(bytes.to_vec())
}
pub fn store(scope: &Scope, secret: &Secret) -> Result<()> {
    let connection = connect()?;
    match item(&connection, scope) {
        Ok(path) => {
            let session = Session::open(connection)?;
            let item = proxy(
                &session.connection,
                path.as_str(),
                "org.freedesktop.Secret.Item",
            )?;
            let data = zeroize::Zeroizing::new(secret.as_bytes().to_vec());
            let _: () = item
                .call_with_flags(
                    "SetSecret",
                    MethodFlags::NoAutoStart.into(),
                    &((
                        session.path.clone(),
                        Vec::<u8>::new(),
                        data.as_slice(),
                        "text/plain; charset=utf-8",
                    ),),
                )
                .map_err(error)?
                .ok_or(Error::Unavailable)?;
            Ok(())
        }
        Err(Error::Missing) => {
            let service = proxy(&connection, ROOT, SERVICE_IFACE)?;
            let collection: OwnedObjectPath = service
                .call_with_flags("ReadAlias", MethodFlags::NoAutoStart.into(), &("default",))
                .map_err(error)?
                .ok_or(Error::Unavailable)?;
            if collection.as_str() == "/" {
                return Err(Error::Unavailable);
            }
            let properties = proxy(
                &connection,
                collection.as_str(),
                "org.freedesktop.DBus.Properties",
            )?;
            let locked: OwnedValue = properties
                .call_with_flags(
                    "Get",
                    MethodFlags::NoAutoStart.into(),
                    &("org.freedesktop.Secret.Collection", "Locked"),
                )
                .map_err(error)?
                .ok_or(Error::Unavailable)?;
            if bool::try_from(locked).map_err(|_| Error::Unavailable)? {
                return Err(Error::Locked);
            }
            drop(properties);
            drop(service);
            let session = Session::open(connection)?;
            let collection = proxy(
                &session.connection,
                collection.as_str(),
                "org.freedesktop.Secret.Collection",
            )?;
            let properties = HashMap::from([
                (
                    "org.freedesktop.Secret.Item.Label",
                    Value::from("Collagene remote access"),
                ),
                (
                    "org.freedesktop.Secret.Item.Attributes",
                    Value::from(attrs(scope)),
                ),
            ]);
            let data = zeroize::Zeroizing::new(secret.as_bytes().to_vec());
            let (created, prompt): (OwnedObjectPath, OwnedObjectPath) = collection
                .call_with_flags(
                    "CreateItem",
                    MethodFlags::NoAutoStart.into(),
                    &(
                        properties,
                        (
                            session.path.clone(),
                            Vec::<u8>::new(),
                            data.as_slice(),
                            "text/plain; charset=utf-8",
                        ),
                        true,
                    ),
                )
                .map_err(error)?
                .ok_or(Error::Unavailable)?;
            if prompt.as_str() != "/" {
                return Err(Error::Locked);
            }
            if created.as_str() == "/" {
                return Err(Error::Unavailable);
            }
            Ok(())
        }
        Err(error) => Err(error),
    }
}
pub fn delete(scope: &Scope) -> Result<()> {
    let connection = connect()?;
    let path = match item(&connection, scope) {
        Ok(path) => path,
        Err(Error::Missing) => return Ok(()),
        Err(error) => return Err(error),
    };
    let item = proxy(&connection, path.as_str(), "org.freedesktop.Secret.Item")?;
    let prompt: OwnedObjectPath = item
        .call_with_flags("Delete", MethodFlags::NoAutoStart.into(), &())
        .map_err(error)?
        .ok_or(Error::Unavailable)?;
    if prompt.as_str() != "/" {
        return Err(Error::Locked);
    }
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn captured_service_owner_must_be_a_unique_bus_peer() {
        assert!(unique_owner(":1.42"));
        assert!(!unique_owner("org.freedesktop.secrets"));
        assert!(!unique_owner(""));
    }
    #[test]
    fn bus_transport_never_executes_or_connects_to_network_fallbacks() {
        for address in ["tcp:host=127.0.0.1,port=1234", "unixexec:path=/bin/false"] {
            assert!(!local_address(&address.parse().unwrap()));
        }
        assert!(local_address(
            &"unix:path=/run/user/1000/bus".parse().unwrap()
        ));
    }
}
