//! Private current/previous archives and independent small preferences. No path
//! comes from IPC; all statements use one fixed host-owned SQLite database.
use crate::{HostError, HostResult};
use collagene_desktop_protocol as wire;
use collagene_library_protocol::{Blob, Digest};
use rusqlite::{Connection, OpenFlags, TransactionBehavior, params};
use std::{path::Path, time::Duration};
const APPLICATION: i64 = 0x43474457;
const SCHEMA: i64 = 1;
const MAX_ARCHIVE: usize = 64 << 20;
pub(crate) struct Working {
    connection: Connection,
}
#[derive(Clone)]
struct Row {
    generation: u64,
    version: Option<Digest>,
    prior: Option<Digest>,
    current: Option<wire::SlotMetadata>,
    previous: Option<wire::SlotMetadata>,
}
fn digest(value: Option<String>) -> HostResult<Option<Digest>> {
    value
        .map(|s| Digest::new(s).map_err(|_| HostError::Corrupt))
        .transpose()
}
fn slot(hash: Option<String>, length: Option<i64>) -> HostResult<Option<wire::SlotMetadata>> {
    match (hash, length) {
        (None, None) => Ok(None),
        (Some(hash), Some(bytes)) if bytes > 0 && bytes <= MAX_ARCHIVE as i64 => {
            Ok(Some(wire::SlotMetadata {
                archive_hash: Digest::new(hash).map_err(|_| HostError::Corrupt)?,
                bytes: bytes as u64,
            }))
        }
        _ => Err(HostError::Corrupt),
    }
}
fn version(
    domain: &str,
    generation: u64,
    prior: &Option<Digest>,
    current: &Option<wire::SlotMetadata>,
    previous: &Option<wire::SlotMetadata>,
) -> HostResult<Digest> {
    let bytes = serde_json::to_vec(&(domain, generation, prior, current, previous))
        .map_err(|_| HostError::Unavailable)?;
    Ok(Digest::of(&bytes))
}
fn row(connection: &Connection) -> HostResult<Row> {
    let raw=connection.query_row("SELECT generation,version,prior_version,current_hash,length(current_archive),previous_hash,length(previous_archive) FROM working WHERE singleton=1",[],|r|Ok((r.get::<_,i64>(0)?,r.get::<_,Option<String>>(1)?,r.get::<_,Option<String>>(2)?,r.get::<_,Option<String>>(3)?,r.get::<_,Option<i64>>(4)?,r.get::<_,Option<String>>(5)?,r.get::<_,Option<i64>>(6)?)))?;
    if raw.0 < 0 {
        return Err(HostError::Corrupt);
    }
    let result = Row {
        generation: raw.0 as u64,
        version: digest(raw.1)?,
        prior: digest(raw.2)?,
        current: slot(raw.3, raw.4)?,
        previous: slot(raw.5, raw.6)?,
    };
    if result.generation == 0 {
        if result.version.is_some()
            || result.prior.is_some()
            || result.current.is_some()
            || result.previous.is_some()
        {
            return Err(HostError::Corrupt);
        }
    } else if result.version.as_ref()
        != Some(&version(
            "desktop-working-1",
            result.generation,
            &result.prior,
            &result.current,
            &result.previous,
        )?)
    {
        return Err(HostError::Corrupt);
    }
    Ok(result)
}
fn metadata(row: &Row) -> wire::WorkingMetadata {
    wire::WorkingMetadata {
        version: row.version.clone(),
        current: row.current.clone(),
        previous: row.previous.clone(),
    }
}
fn check(row: &Row, expected: &Option<Digest>) -> HostResult<()> {
    if &row.version != expected {
        return Err(HostError::WorkingConflict(metadata(row)));
    }
    Ok(())
}
fn validate(bytes: &[u8]) -> HostResult<wire::SlotMetadata> {
    if bytes.is_empty() || bytes.len() > MAX_ARCHIVE {
        return Err(HostError::Invalid(
            "Native workspace archive exceeds its slot bound.",
        ));
    }
    collagene_core::workspace::load_workspace(bytes)
        .map_err(|_| HostError::Invalid("Invalid native workspace archive."))?;
    Ok(wire::SlotMetadata {
        archive_hash: Digest::of(bytes),
        bytes: bytes.len() as u64,
    })
}
fn read_slot(
    connection: &Connection,
    slot: wire::WorkingSlot,
    expected: &Option<Digest>,
) -> HostResult<(Digest, Vec<u8>)> {
    let info = row(connection)?;
    check(&info, expected)?;
    let descriptor = match slot {
        wire::WorkingSlot::Current => info.current.as_ref(),
        wire::WorkingSlot::Previous => info.previous.as_ref(),
    }
    .ok_or(HostError::Invalid("The requested working slot is empty."))?;
    let sql = match slot {
        wire::WorkingSlot::Current => {
            "SELECT CASE WHEN length(current_archive)<=67108864 THEN current_archive ELSE NULL END FROM working WHERE singleton=1"
        }
        wire::WorkingSlot::Previous => {
            "SELECT CASE WHEN length(previous_archive)<=67108864 THEN previous_archive ELSE NULL END FROM working WHERE singleton=1"
        }
    };
    let bytes: Option<Vec<u8>> = connection.query_row(sql, [], |r| r.get(0))?;
    let bytes = bytes.ok_or(HostError::Corrupt)?;
    if Digest::of(&bytes) != descriptor.archive_hash || bytes.len() as u64 != descriptor.bytes {
        return Err(HostError::Corrupt);
    }
    validate(&bytes).map_err(|_| HostError::Corrupt)?;
    Ok((info.version.ok_or(HostError::Corrupt)?, bytes))
}
impl Working {
    pub(crate) fn open(path: &Path) -> HostResult<Self> {
        let connection = Connection::open_with_flags(
            path,
            OpenFlags::SQLITE_OPEN_READ_WRITE
                | OpenFlags::SQLITE_OPEN_CREATE
                | OpenFlags::SQLITE_OPEN_NOFOLLOW,
        )?;
        connection.busy_timeout(Duration::from_secs(5))?;
        let application: i64 = connection.query_row("PRAGMA application_id", [], |r| r.get(0))?;
        let schema: i64 = connection.query_row("PRAGMA user_version", [], |r| r.get(0))?;
        if application == 0 && schema == 0 {
            let tables: i64 = connection.query_row(
                "SELECT count(*) FROM sqlite_master WHERE name NOT LIKE 'sqlite_%'",
                [],
                |r| r.get(0),
            )?;
            if tables != 0 {
                return Err(HostError::Corrupt);
            }
        } else if application != APPLICATION || schema != SCHEMA {
            return Err(HostError::Corrupt);
        }
        connection.execute_batch("PRAGMA journal_mode=WAL; PRAGMA synchronous=FULL; PRAGMA foreign_keys=ON;
            BEGIN IMMEDIATE;
            CREATE TABLE IF NOT EXISTS working(singleton INTEGER PRIMARY KEY CHECK(singleton=1),generation INTEGER NOT NULL CHECK(generation>=0),version TEXT,prior_version TEXT,current_hash TEXT,current_archive BLOB,previous_hash TEXT,previous_archive BLOB,CHECK((current_hash IS NULL)=(current_archive IS NULL)),CHECK((previous_hash IS NULL)=(previous_archive IS NULL)),CHECK(current_archive IS NULL OR length(current_archive)<=67108864),CHECK(previous_archive IS NULL OR length(previous_archive)<=67108864));
            INSERT OR IGNORE INTO working VALUES(1,0,NULL,NULL,NULL,NULL,NULL,NULL);
            CREATE TABLE IF NOT EXISTS preferences(singleton INTEGER PRIMARY KEY CHECK(singleton=1),generation INTEGER NOT NULL CHECK(generation>=0),version TEXT,prior_version TEXT,bytes BLOB NOT NULL CHECK(length(bytes)<=1024));
            INSERT OR IGNORE INTO preferences VALUES(1,0,NULL,NULL,x'');
            PRAGMA application_id=1128744023; PRAGMA user_version=1; COMMIT;")?;
        let result = Self { connection };
        result.metadata()?;
        result.preferences()?;
        Ok(result)
    }
    pub(crate) fn metadata(&self) -> HostResult<wire::WorkingMetadata> {
        Ok(metadata(&row(&self.connection)?))
    }
    pub(crate) fn read(
        &mut self,
        slot: wire::WorkingSlot,
        expected: &Option<Digest>,
    ) -> HostResult<(Digest, Vec<u8>)> {
        let tx = self.connection.transaction()?;
        let result = read_slot(&tx, slot, expected)?;
        tx.commit()?;
        Ok(result)
    }
    pub(crate) fn save<F>(
        &mut self,
        expected: &Option<Digest>,
        current: &[u8],
        previous: Option<&[u8]>,
        gate: F,
    ) -> HostResult<wire::WorkingMetadata>
    where
        F: FnOnce() -> HostResult<()>,
    {
        let current_meta = validate(current)?;
        let previous_meta = previous.map(validate).transpose()?;
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        let old = row(&tx)?;
        check(&old, expected)?;
        let previous_meta = if previous.is_some() {
            previous_meta
        } else {
            old.previous
        };
        let generation = old
            .generation
            .checked_add(1)
            .filter(|n| *n <= i64::MAX as u64)
            .ok_or(HostError::Unavailable)?;
        let version = version(
            "desktop-working-1",
            generation,
            &old.version,
            &Some(current_meta.clone()),
            &previous_meta,
        )?;
        gate()?;
        if let Some(previous) = previous {
            tx.execute("UPDATE working SET generation=?1,version=?2,prior_version=?3,current_hash=?4,current_archive=?5,previous_hash=?6,previous_archive=?7 WHERE singleton=1",params![generation as i64,version.as_str(),old.version.as_ref().map(Digest::as_str),current_meta.archive_hash.as_str(),current,previous_meta.as_ref().map(|m|m.archive_hash.as_str()),previous])?;
        } else {
            tx.execute("UPDATE working SET generation=?1,version=?2,prior_version=?3,current_hash=?4,current_archive=?5 WHERE singleton=1",params![generation as i64,version.as_str(),old.version.as_ref().map(Digest::as_str),current_meta.archive_hash.as_str(),current])?;
        }
        let result = metadata(&row(&tx)?);
        tx.commit()?;
        Ok(result)
    }
    pub(crate) fn recover<F>(
        &mut self,
        expected: &Digest,
        gate: F,
    ) -> HostResult<wire::WorkingMetadata>
    where
        F: FnOnce() -> HostResult<()>,
    {
        let (_, previous) = self.read(wire::WorkingSlot::Previous, &Some(expected.clone()))?;
        self.save(&Some(expected.clone()), &previous, None, gate)
    }
    pub(crate) fn preferences(&self) -> HostResult<wire::Preferences> {
        let(generation,version_raw,prior_raw,bytes):(i64,Option<String>,Option<String>,Option<Vec<u8>>)=self.connection.query_row("SELECT generation,version,prior_version,CASE WHEN length(bytes)<=1024 THEN bytes ELSE NULL END FROM preferences WHERE singleton=1",[],|r|Ok((r.get(0)?,r.get(1)?,r.get(2)?,r.get(3)?)))?;
        let bytes = bytes.ok_or(HostError::Corrupt)?;
        let version_value = digest(version_raw)?;
        let prior = digest(prior_raw)?;
        if generation == 0 {
            if version_value.is_some() || prior.is_some() || !bytes.is_empty() {
                return Err(HostError::Corrupt);
            }
        } else if generation < 0
            || version_value.as_ref()
                != Some(&Digest::of(
                    &serde_json::to_vec(&(
                        "desktop-preferences-1",
                        generation,
                        prior,
                        Digest::of(&bytes),
                    ))
                    .map_err(|_| HostError::Unavailable)?,
                ))
        {
            return Err(HostError::Corrupt);
        }
        Ok(wire::Preferences {
            version: version_value,
            bytes: Blob::new(bytes).map_err(|_| HostError::Corrupt)?,
        })
    }
    pub(crate) fn save_preferences<F>(
        &mut self,
        expected: &Option<Digest>,
        bytes: &[u8],
        gate: F,
    ) -> HostResult<wire::Preferences>
    where
        F: FnOnce() -> HostResult<()>,
    {
        if bytes.len() > 1024 {
            return Err(HostError::Invalid("View preferences exceed1024 bytes."));
        }
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        let (generation, raw): (i64, Option<String>) = tx.query_row(
            "SELECT generation,version FROM preferences WHERE singleton=1",
            [],
            |r| Ok((r.get(0)?, r.get(1)?)),
        )?;
        let prior = digest(raw)?;
        if &prior != expected {
            drop(tx);
            return Err(HostError::PreferencesConflict(self.preferences()?));
        }
        let generation = generation.checked_add(1).ok_or(HostError::Unavailable)?;
        let new = Digest::of(
            &serde_json::to_vec(&(
                "desktop-preferences-1",
                generation,
                &prior,
                Digest::of(bytes),
            ))
            .map_err(|_| HostError::Unavailable)?,
        );
        gate()?;
        tx.execute("UPDATE preferences SET generation=?1,version=?2,prior_version=?3,bytes=?4 WHERE singleton=1",params![generation,new.as_str(),prior.as_ref().map(Digest::as_str),bytes])?;
        tx.commit()?;
        Ok(wire::Preferences {
            version: Some(new),
            bytes: Blob::new(bytes.to_vec()).map_err(|_| HostError::Unavailable)?,
        })
    }
}
