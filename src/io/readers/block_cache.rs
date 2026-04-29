use std::sync::Arc;

use quick_cache::sync::Cache;

pub(crate) const DEFAULT_TOTAL_BLOCK_CACHE_ENTRIES: usize = 500;

type BlockCacheKey = (u64, u64);

/// Split the total block-cache budget across `sample_count` files.
///
/// The first `total_capacity % sample_count` files each get one extra entry
/// so that the per-file capacities sum to exactly `total_capacity`.
pub fn split_block_cache_capacity(
    total_capacity: usize,
    sample_count: usize,
    sample_index: usize,
) -> usize {
    if sample_count == 0 || sample_index >= sample_count {
        return 0;
    }
    let base = total_capacity / sample_count;
    let remainder = total_capacity % sample_count;
    if sample_index < remainder {
        base + 1
    } else {
        base
    }
}

/// A thread-safe shared block cache with LRU eviction.
///
/// Multiple worker threads share a single instance via `Arc<SharedBlockCache>`.
/// When the cache reaches capacity, the least recently used entry is evicted.
/// A zero-capacity cache acts as a no-op (never stores entries).
pub struct SharedBlockCache {
    cache: Option<Cache<BlockCacheKey, Arc<[u8]>>>,
}

impl SharedBlockCache {
    pub fn new() -> Self {
        Self::with_capacity(DEFAULT_TOTAL_BLOCK_CACHE_ENTRIES)
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            cache: (capacity > 0).then(|| Cache::new(capacity)),
        }
    }

    pub fn get(&self, key: &BlockCacheKey) -> Option<Arc<[u8]>> {
        self.cache.as_ref().and_then(|cache| cache.get(key))
    }

    pub fn insert(&self, key: BlockCacheKey, value: Arc<[u8]>) {
        if let Some(cache) = &self.cache {
            cache.insert(key, value);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn split_cache_capacity_divides_total_with_remainder() {
        let caps: Vec<_> = (0..3)
            .map(|sample_index| split_block_cache_capacity(500, 3, sample_index))
            .collect();
        assert_eq!(caps, vec![167, 167, 166]);
        assert_eq!(caps.iter().sum::<usize>(), 500);
    }

    #[test]
    fn split_cache_capacity_allows_zero_capacity_when_samples_exceed_total() {
        assert_eq!(split_block_cache_capacity(3, 5, 0), 1);
        assert_eq!(split_block_cache_capacity(3, 5, 1), 1);
        assert_eq!(split_block_cache_capacity(3, 5, 2), 1);
        assert_eq!(split_block_cache_capacity(3, 5, 3), 0);
        assert_eq!(split_block_cache_capacity(3, 5, 4), 0);
    }

    #[test]
    fn zero_capacity_cache_never_stores_entries() {
        let cache = SharedBlockCache::with_capacity(0);
        let key = (128, 64);
        cache.insert(key, Arc::from(&b"payload"[..]));
        assert!(cache.get(&key).is_none());
    }

    #[test]
    fn independent_per_file_caches_do_not_share_offsets() {
        let cache_a = SharedBlockCache::with_capacity(2);
        let cache_b = SharedBlockCache::with_capacity(2);
        let key = (128, 64);

        cache_a.insert(key, Arc::from(&b"file-a"[..]));
        cache_b.insert(key, Arc::from(&b"file-b"[..]));

        assert_eq!(cache_a.get(&key).unwrap().as_ref(), b"file-a");
        assert_eq!(cache_b.get(&key).unwrap().as_ref(), b"file-b");
    }
}
