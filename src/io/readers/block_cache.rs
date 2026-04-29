use std::sync::Arc;

use quick_cache::sync::Cache;

pub(crate) const DEFAULT_PER_FILE_BLOCK_CACHE_ENTRIES: usize = 200;
pub(crate) const HARD_LIMIT_TOTAL_BLOCK_CACHE_ENTRIES: usize = 2000;

type BlockCacheKey = (u64, u64);

/// Compute the block-cache capacity for file `sample_index` out of `file_count`.
///
/// Each file gets up to `DEFAULT_PER_FILE_BLOCK_CACHE_ENTRIES` (200) entries.
/// If `file_count * 200` exceeds the hard limit (2000), the 2000-entry budget
/// is distributed with remainder: the first `remainder` files each get
/// `base + 1`, the rest get `base` (where `base = 2000 / file_count`).
/// This ensures `sum(per_file) == min(file_count * 200, 2000)` exactly.
pub fn compute_per_file_block_cache_capacity(file_count: usize, sample_index: usize) -> usize {
    if file_count == 0 || sample_index >= file_count {
        return 0;
    }
    let desired = file_count.saturating_mul(DEFAULT_PER_FILE_BLOCK_CACHE_ENTRIES);
    if desired <= HARD_LIMIT_TOTAL_BLOCK_CACHE_ENTRIES {
        DEFAULT_PER_FILE_BLOCK_CACHE_ENTRIES
    } else {
        let base = HARD_LIMIT_TOTAL_BLOCK_CACHE_ENTRIES / file_count;
        let remainder = HARD_LIMIT_TOTAL_BLOCK_CACHE_ENTRIES % file_count;
        if sample_index < remainder {
            base + 1
        } else {
            base
        }
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
        Self::with_capacity(DEFAULT_PER_FILE_BLOCK_CACHE_ENTRIES)
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
    fn few_files_get_default_per_file() {
        for i in 0..3 {
            assert_eq!(compute_per_file_block_cache_capacity(3, i), 200);
        }
    }

    #[test]
    fn exactly_at_hard_limit() {
        for i in 0..10 {
            assert_eq!(compute_per_file_block_cache_capacity(10, i), 200);
        }
    }

    #[test]
    fn over_hard_limit_distributes_with_remainder() {
        // 11 files: budget=2000, 2000/11=181 base, remainder=9
        // files 0..8 get 182, files 9..10 get 181
        assert_eq!(compute_per_file_block_cache_capacity(11, 0), 182);
        assert_eq!(compute_per_file_block_cache_capacity(11, 8), 182);
        assert_eq!(compute_per_file_block_cache_capacity(11, 9), 181);
        assert_eq!(compute_per_file_block_cache_capacity(11, 10), 181);
        // 20 files: 2000/20=100, remainder=0, all get 100
        for i in 0..20 {
            assert_eq!(compute_per_file_block_cache_capacity(20, i), 100);
        }
    }

    #[test]
    fn many_files_first_2000_get_one_rest_zero() {
        // 3000 files: 2000/3000=0 base, remainder=2000
        // first 2000 files get 1, last 1000 get 0
        assert_eq!(compute_per_file_block_cache_capacity(3000, 0), 1);
        assert_eq!(compute_per_file_block_cache_capacity(3000, 1999), 1);
        assert_eq!(compute_per_file_block_cache_capacity(3000, 2000), 0);
        assert_eq!(compute_per_file_block_cache_capacity(3000, 2999), 0);
    }

    #[test]
    fn total_equals_budget_exactly() {
        for file_count in [1, 3, 10, 11, 20, 50, 100, 500, 2000, 3000] {
            let total: usize = (0..file_count)
                .map(|i| compute_per_file_block_cache_capacity(file_count, i))
                .sum();
            let expected = std::cmp::min(
                file_count * DEFAULT_PER_FILE_BLOCK_CACHE_ENTRIES,
                HARD_LIMIT_TOTAL_BLOCK_CACHE_ENTRIES,
            );
            assert_eq!(
                total, expected,
                "file_count={file_count}: total={total} != expected={expected}"
            );
        }
    }

    #[test]
    fn single_file_gets_default() {
        assert_eq!(compute_per_file_block_cache_capacity(1, 0), 200);
    }

    #[test]
    fn zero_files_returns_zero() {
        assert_eq!(compute_per_file_block_cache_capacity(0, 0), 0);
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
