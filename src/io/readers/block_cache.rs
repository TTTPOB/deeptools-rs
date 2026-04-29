use std::sync::Arc;

use quick_cache::sync::Cache;

const MAX_SHARED_BLOCK_CACHE_ENTRIES: usize = 500;

/// A thread-safe shared block cache with LRU eviction.
///
/// Multiple worker threads share a single instance via `Arc<SharedBlockCache>`.
/// When the cache reaches capacity, the least recently used entry is evicted.
pub struct SharedBlockCache {
    cache: Cache<(u64, u64), Arc<[u8]>>,
}

impl SharedBlockCache {
    pub fn new() -> Self {
        Self {
            cache: Cache::new(MAX_SHARED_BLOCK_CACHE_ENTRIES),
        }
    }

    pub fn get(&self, key: &(u64, u64)) -> Option<Arc<[u8]>> {
        self.cache.get(key)
    }

    pub fn insert(&self, key: (u64, u64), value: Arc<[u8]>) {
        self.cache.insert(key, value);
    }
}
