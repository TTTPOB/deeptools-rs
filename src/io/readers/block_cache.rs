use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use dashmap::DashMap;

const MAX_SHARED_BLOCK_CACHE_ENTRIES: usize = 500;

/// A thread-safe shared block cache backed by DashMap.
///
/// Multiple worker threads share a single instance via `Arc<SharedBlockCache>`.
/// When the cache reaches its capacity limit, new inserts are silently skipped
/// ("skip insert when full" policy) rather than evicting existing entries.
pub struct SharedBlockCache {
    map: DashMap<(u64, u64), Arc<[u8]>>,
    len: AtomicUsize,
}

impl SharedBlockCache {
    pub fn new() -> Self {
        Self {
            map: DashMap::with_capacity(MAX_SHARED_BLOCK_CACHE_ENTRIES),
            len: AtomicUsize::new(0),
        }
    }

    pub fn get(&self, key: &(u64, u64)) -> Option<Arc<[u8]>> {
        self.map.get(key).map(|entry| Arc::clone(entry.value()))
    }

    /// Insert if not full and not already present. Skip silently if full.
    ///
    /// Note: the `AtomicUsize` len counter may drift slightly under high
    /// concurrency (two threads both check len < 500, both insert). This is
    /// acceptable — we allow slightly exceeding the cap.
    pub fn insert(&self, key: (u64, u64), value: Arc<[u8]>) {
        if self.len.load(Ordering::Relaxed) >= MAX_SHARED_BLOCK_CACHE_ENTRIES {
            return;
        }
        use dashmap::mapref::entry::Entry;
        match self.map.entry(key) {
            Entry::Occupied(_) => {}
            Entry::Vacant(v) => {
                v.insert(value);
                self.len.fetch_add(1, Ordering::Relaxed);
            }
        }
    }
}
