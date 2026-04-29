use serde_json::Value;

/// A single header key/value difference.
#[derive(Debug, Clone)]
pub struct HeaderDiff {
    pub key: String,
    pub left: String,
    pub right: String,
}

impl HeaderDiff {
    pub fn new(key: impl Into<String>, left: &Value, right: &Value) -> Self {
        Self {
            key: key.into(),
            left: value_display(left),
            right: value_display(right),
        }
    }
}

/// Compare two JSON header objects key-by-key.
///
/// Keys listed in `ignore_keys` are skipped entirely.
/// Only keys present in either object are compared; a missing key in one side
/// is reported as a difference (shown as `<missing>`).
pub fn compare_headers(left: &Value, right: &Value, ignore_keys: &[String]) -> Vec<HeaderDiff> {
    let mut diffs = Vec::new();

    let left_obj = match left.as_object() {
        Some(o) => o,
        None => {
            // Treat the entire header as a single diff if not an object
            if left != right {
                diffs.push(HeaderDiff {
                    key: "<root>".into(),
                    left: value_display(left),
                    right: value_display(right),
                });
            }
            return diffs;
        }
    };
    let right_obj = match right.as_object() {
        Some(o) => o,
        None => {
            diffs.push(HeaderDiff {
                key: "<root>".into(),
                left: value_display(left),
                right: value_display(right),
            });
            return diffs;
        }
    };

    // Collect all keys from both sides, in stable order (left first, then right-only)
    let mut all_keys: Vec<&str> = left_obj.keys().map(|k| k.as_str()).collect();
    for key in right_obj.keys() {
        if !left_obj.contains_key(key.as_str()) {
            all_keys.push(key.as_str());
        }
    }

    for key in all_keys {
        if ignore_keys.iter().any(|k| k == key) {
            continue;
        }

        match (left_obj.get(key), right_obj.get(key)) {
            (Some(l), Some(r)) => {
                if l != r {
                    diffs.push(HeaderDiff::new(key, l, r));
                }
            }
            (Some(l), None) => {
                diffs.push(HeaderDiff {
                    key: key.to_owned(),
                    left: value_display(l),
                    right: "<missing>".into(),
                });
            }
            (None, Some(r)) => {
                diffs.push(HeaderDiff {
                    key: key.to_owned(),
                    left: "<missing>".into(),
                    right: value_display(r),
                });
            }
            (None, None) => {} // impossible — key came from one of them
        }
    }

    diffs
}

fn value_display(v: &Value) -> String {
    match v {
        Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}
