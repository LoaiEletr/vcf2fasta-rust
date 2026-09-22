//! Pinned (page-locked) host memory for GPU staging.

use std::marker::PhantomData;
use std::ops::{Deref, DerefMut};

enum Storage<T> {
    // Constructed only when the `cuda` feature is enabled; allow the
    // variant to exist unconditionally so the type signature is stable.
    #[allow(dead_code)]
    Pinned { ptr: *mut T, len: usize },
    Regular(Vec<T>),
}

pub struct PinnedBuf<T: Copy + Default> {
    storage: Storage<T>,
    _marker: PhantomData<T>,
}

impl<T: Copy + Default> PinnedBuf<T> {
    pub fn new(len: usize) -> Self {
        #[cfg(feature = "cuda")]
        {
            if len > 0 {
                let elem_size = std::mem::size_of::<T>().max(1);
                if let Some(bytes) = len.checked_mul(elem_size) {
                    if bytes > 0 {
                        let raw = unsafe {
                            cudarc::driver::result::malloc_host(bytes, 0u32)
                        };
                        if let Ok(raw) = raw {
                            if !raw.is_null() {
                                let ptr = raw as *mut T;
                                unsafe {
                                    for i in 0..len {
                                        std::ptr::write(ptr.add(i), T::default());
                                    }
                                }
                                return Self {
                                    storage: Storage::Pinned { ptr, len },
                                    _marker: PhantomData,
                                };
                            }
                        }
                    }
                }
            }
        }
        Self {
            storage: Storage::Regular(vec![T::default(); len]),
            _marker: PhantomData,
        }
    }

    /// Grows the buffer to at least `len` elements.
    ///
    /// Reallocates only when growing. Never shrinks. This is the whole
    /// point of `BatchScratch`: amortise `cudaHostAlloc` across many
    /// batches by only reallocating when a batch is bigger than any
    /// seen before.
    pub fn ensure_capacity(&mut self, len: usize) {
        if len > self.len() {
            *self = Self::new(len);
        }
    }

    pub fn is_pinned(&self) -> bool {
        matches!(self.storage, Storage::Pinned { .. })
    }

    pub fn len(&self) -> usize {
        match &self.storage {
            Storage::Pinned { len, .. } => *len,
            Storage::Regular(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<T: Copy + Default> Deref for PinnedBuf<T> {
    type Target = [T];
    fn deref(&self) -> &[T] {
        match &self.storage {
            Storage::Pinned { ptr, len } => unsafe {
                std::slice::from_raw_parts(*ptr, *len)
            },
            Storage::Regular(v) => v.as_slice(),
        }
    }
}

impl<T: Copy + Default> DerefMut for PinnedBuf<T> {
    fn deref_mut(&mut self) -> &mut [T] {
        match &mut self.storage {
            Storage::Pinned { ptr, len } => unsafe {
                std::slice::from_raw_parts_mut(*ptr, *len)
            },
            Storage::Regular(v) => v.as_mut_slice(),
        }
    }
}

impl<T: Copy + Default> Drop for PinnedBuf<T> {
    fn drop(&mut self) {
        if let Storage::Pinned { ptr, .. } = self.storage {
            #[cfg(feature = "cuda")]
            unsafe {
                let _ = cudarc::driver::result::free_host(
                    ptr as *mut std::ffi::c_void,
                );
            }
            #[cfg(not(feature = "cuda"))]
            {
                let _ = ptr;
            }
        }
    }
}

unsafe impl<T: Copy + Default + Send> Send for PinnedBuf<T> {}
unsafe impl<T: Copy + Default + Sync> Sync for PinnedBuf<T> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pinned_buf_allocates_and_zeroes() {
        let b: PinnedBuf<u8> = PinnedBuf::new(1024);
        assert_eq!(b.len(), 1024);
        assert!(b.iter().all(|&x| x == 0));
    }

    #[test]
    fn ensure_capacity_grows() {
        let mut b: PinnedBuf<u8> = PinnedBuf::new(16);
        b.ensure_capacity(1024);
        assert!(b.len() >= 1024);
    }

    #[test]
    fn ensure_capacity_never_shrinks() {
        let mut b: PinnedBuf<u8> = PinnedBuf::new(1024);
        b.ensure_capacity(16);
        assert_eq!(b.len(), 1024);
    }

    #[test]
    fn ensure_capacity_is_noop_when_big_enough() {
        let mut b: PinnedBuf<u8> = PinnedBuf::new(1024);
        let ptr_before = b.as_ptr();
        b.ensure_capacity(512);
        assert_eq!(b.as_ptr(), ptr_before);
    }
}