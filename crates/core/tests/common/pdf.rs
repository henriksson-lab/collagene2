use std::collections::BTreeMap;
/// Independent narrow PDF reader for these uncompressed development outputs.
/// It decodes literal strings and the PDF ToUnicode mapping, not exporter internals.
fn literals(bytes: &[u8]) -> Vec<(Vec<u8>, usize)> {
    let mut result = vec![];
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] != b'(' {
            i += 1;
            continue;
        }
        let mut out = vec![];
        let mut depth = 1;
        i += 1;
        while i < bytes.len() && depth > 0 {
            let b = bytes[i];
            i += 1;
            match b {
                b'\\' => {
                    if i >= bytes.len() {
                        break;
                    }
                    let escaped = bytes[i];
                    i += 1;
                    if escaped.is_ascii_digit() && escaped <= b'7' {
                        let mut n = escaped - b'0';
                        for _ in 0..2 {
                            if i < bytes.len() && (b'0'..=b'7').contains(&bytes[i]) {
                                n = n.wrapping_mul(8).wrapping_add(bytes[i] - b'0');
                                i += 1;
                            } else {
                                break;
                            }
                        }
                        out.push(n);
                    } else {
                        match escaped {
                            b'n' => out.push(b'\n'),
                            b'r' => out.push(b'\r'),
                            b't' => out.push(b'\t'),
                            b'b' => out.push(8),
                            b'f' => out.push(12),
                            b'\n' => {}
                            b'\r' => {
                                if bytes.get(i) == Some(&b'\n') {
                                    i += 1;
                                }
                            }
                            x => out.push(x),
                        }
                    }
                }
                b'(' => {
                    depth += 1;
                    out.push(b);
                }
                b')' => {
                    depth -= 1;
                    if depth > 0 {
                        out.push(b);
                    }
                }
                _ => out.push(b),
            }
        }
        if depth == 0 {
            result.push((out, i));
        }
    }
    result
}
pub fn extracted(bytes: &[u8]) -> String {
    let text = String::from_utf8_lossy(bytes);
    let mut map = BTreeMap::<u16, String>::new();
    for line in text.lines() {
        let tokens: Vec<_> = line.split_whitespace().collect();
        if tokens.len() != 2
            || !tokens[0].starts_with('<')
            || tokens[0].len() != 6
            || !tokens[1].starts_with('<')
        {
            continue;
        }
        let Ok(cid) = u16::from_str_radix(&tokens[0][1..5], 16) else {
            continue;
        };
        let unicode = &tokens[1][1..tokens[1].len() - 1];
        if unicode.len() % 4 != 0 {
            continue;
        }
        let units: Option<Vec<_>> = (0..unicode.len())
            .step_by(4)
            .map(|i| u16::from_str_radix(&unicode[i..i + 4], 16).ok())
            .collect();
        if let Some(units) = units {
            if let Ok(s) = String::from_utf16(&units) {
                map.insert(cid, s);
            }
        }
    }
    let mut out = String::new();
    // Delimit real streams by their declared /Length, independent of embedded bytes.
    let mut start = 0;
    while let Some(relative) = bytes[start..].windows(8).position(|s| s == b"\nstream\n") {
        let content = start + relative + 8;
        let dictionary = String::from_utf8_lossy(&bytes[start..content]);
        let Some(length) = dictionary
            .rsplit("/Length ")
            .next()
            .and_then(|s| s.split_whitespace().next())
            .and_then(|s| s.parse::<usize>().ok())
        else {
            break;
        };
        if content + length > bytes.len() {
            break;
        }
        let stream = &bytes[content..content + length];
        if stream.windows(3).any(|s| s == b"BT\n") {
            for (s, end) in literals(stream) {
                if !stream[end..]
                    .iter()
                    .skip_while(|b| b.is_ascii_whitespace())
                    .copied()
                    .take(2)
                    .eq(b"Tj".iter().copied())
                {
                    continue;
                }
                for pair in s.chunks_exact(2) {
                    if let Some(t) = map.get(&u16::from_be_bytes([pair[0], pair[1]])) {
                        out.push_str(t);
                    }
                }
                out.push('\n');
            }
        }
        start = content + length;
    }
    out
}
pub fn check_xref(bytes: &[u8]) {
    let text = String::from_utf8_lossy(bytes);
    let start = text
        .rsplit("startxref\n")
        .next()
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .parse::<usize>()
        .unwrap();
    assert_eq!(&bytes[start..start + 5], b"xref\n");
    let xref = std::str::from_utf8(&bytes[start..]).unwrap();
    let mut lines = xref.lines();
    assert_eq!(lines.next(), Some("xref"));
    let header = lines.next().unwrap().split_whitespace().collect::<Vec<_>>();
    assert_eq!(header[0], "0");
    let count = header[1].parse::<usize>().unwrap();
    for id in 0..count {
        let row = lines.next().unwrap();
        if id == 0 {
            assert_eq!(row.split_whitespace().last(), Some("f"));
            continue;
        }
        let offset = row[..10].parse::<usize>().unwrap();
        assert!(
            bytes[offset..].starts_with(format!("{id} 0 obj\n").as_bytes()),
            "xref {id} -> {offset}"
        );
    }
    assert!(bytes.ends_with(b"%%EOF") || bytes.ends_with(b"%%EOF\n"));
}
