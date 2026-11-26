//! Deserialize borrowed binary lengths before allocating owned storage.
use serde::{
    Deserialize, Deserializer, Serialize, Serializer,
    de::{Error, SeqAccess, Visitor},
};
use std::{fmt, ops::Deref};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Text<const N: usize>(String);
impl<const N: usize> Text<N> {
    pub fn new(value: impl AsRef<str>) -> Result<Self, crate::Error> {
        let s = value.as_ref();
        if s.len() > N {
            return Err(crate::Error::TooLarge);
        }
        Ok(Self(s.to_owned()))
    }
    pub fn as_str(&self) -> &str {
        &self.0
    }
    pub fn into_string(self) -> String {
        self.0
    }
}
impl<const N: usize> Deref for Text<N> {
    type Target = str;
    fn deref(&self) -> &str {
        &self.0
    }
}
impl<const N: usize> Serialize for Text<N> {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(&self.0)
    }
}
impl<'de, const N: usize> Deserialize<'de> for Text<N> {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V<const N: usize>;
        impl<'de, const N: usize> Visitor<'de> for V<N> {
            type Value = Text<N>;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "at most {N} UTF-8 bytes")
            }
            fn visit_str<E: Error>(self, s: &str) -> Result<Self::Value, E> {
                if s.len() > N {
                    return Err(E::custom("text bound exceeded"));
                }
                Ok(Text(s.to_owned()))
            }
            fn visit_borrowed_str<E: Error>(self, s: &'de str) -> Result<Self::Value, E> {
                self.visit_str(s)
            }
            fn visit_string<E: Error>(self, s: String) -> Result<Self::Value, E> {
                if s.len() > N {
                    return Err(E::custom("text bound exceeded"));
                }
                Ok(Text(s))
            }
        }
        d.deserialize_str(V::<N>)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Blob<const N: usize>(Vec<u8>);
impl<const N: usize> Blob<N> {
    pub fn new(bytes: Vec<u8>) -> Result<Self, crate::Error> {
        if bytes.len() > N {
            return Err(crate::Error::TooLarge);
        }
        Ok(Self(bytes))
    }
    pub fn as_slice(&self) -> &[u8] {
        &self.0
    }
    pub fn into_vec(self) -> Vec<u8> {
        self.0
    }
}
impl<const N: usize> Deref for Blob<N> {
    type Target = [u8];
    fn deref(&self) -> &[u8] {
        &self.0
    }
}
impl<const N: usize> Serialize for Blob<N> {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_bytes(&self.0)
    }
}
impl<'de, const N: usize> Deserialize<'de> for Blob<N> {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V<const N: usize>;
        impl<'de, const N: usize> Visitor<'de> for V<N> {
            type Value = Blob<N>;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "at most {N} binary bytes")
            }
            fn visit_bytes<E: Error>(self, b: &[u8]) -> Result<Self::Value, E> {
                if b.len() > N {
                    return Err(E::custom("blob bound exceeded"));
                }
                Ok(Blob(b.to_vec()))
            }
            fn visit_borrowed_bytes<E: Error>(self, b: &'de [u8]) -> Result<Self::Value, E> {
                self.visit_bytes(b)
            }
            fn visit_byte_buf<E: Error>(self, b: Vec<u8>) -> Result<Self::Value, E> {
                if b.len() > N {
                    return Err(E::custom("blob bound exceeded"));
                }
                Ok(Blob(b))
            }
            fn visit_seq<A: SeqAccess<'de>>(self, mut a: A) -> Result<Self::Value, A::Error> {
                if a.size_hint().is_some_and(|n| n > N) {
                    return Err(A::Error::custom("blob bound exceeded"));
                }
                let mut b = Vec::with_capacity(a.size_hint().unwrap_or(0).min(N));
                while let Some(x) = a.next_element()? {
                    if b.len() == N {
                        return Err(A::Error::custom("blob bound exceeded"));
                    }
                    b.push(x);
                }
                Ok(Blob(b))
            }
        }
        d.deserialize_bytes(V::<N>)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Rows<T, const N: usize>(Vec<T>);
impl<T, const N: usize> Rows<T, N> {
    pub fn new(rows: Vec<T>) -> Result<Self, crate::Error> {
        if rows.len() > N {
            return Err(crate::Error::TooLarge);
        }
        Ok(Self(rows))
    }
    pub fn as_slice(&self) -> &[T] {
        &self.0
    }
    pub fn into_vec(self) -> Vec<T> {
        self.0
    }
}
impl<T, const N: usize> Deref for Rows<T, N> {
    type Target = [T];
    fn deref(&self) -> &[T] {
        &self.0
    }
}
impl<T: Serialize, const N: usize> Serialize for Rows<T, N> {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        self.0.serialize(s)
    }
}
impl<'de, T: Deserialize<'de>, const N: usize> Deserialize<'de> for Rows<T, N> {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V<T, const N: usize>(std::marker::PhantomData<T>);
        impl<'de, T: Deserialize<'de>, const N: usize> Visitor<'de> for V<T, N> {
            type Value = Rows<T, N>;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "at most {N} rows")
            }
            fn visit_seq<A: SeqAccess<'de>>(self, mut a: A) -> Result<Self::Value, A::Error> {
                if a.size_hint().is_some_and(|n| n > N) {
                    return Err(A::Error::custom("row bound exceeded"));
                }
                let mut rows = Vec::with_capacity(a.size_hint().unwrap_or(0).min(N));
                while let Some(row) = a.next_element()? {
                    if rows.len() == N {
                        return Err(A::Error::custom("row bound exceeded"));
                    }
                    rows.push(row);
                }
                Ok(Rows(rows))
            }
        }
        d.deserialize_seq(V::<T, N>(std::marker::PhantomData))
    }
}
