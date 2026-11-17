module atmos_phys_kiss_rng
   ! The KISS (Keep It Simple Stupid) random number generator for one stream,
   ! callable from OpenACC kernels.
   !
   ! kiss_rng_next gives the same numbers as one stream of the share code's
   ! kissvec (share/RandNum/src/kissvec/kissvec.c, used by ShrKissRandGen):
   ! the four seeds advance with the same unsigned 32-bit arithmetic, emulated
   ! in 64-bit integers masked to 32 bits (Fortran has no unsigned type, and
   ! signed overflow is not portable), and the result is scaled as there:
   ! int32(sum of seeds) * 2.328306e-10 + 0.5. Each stream (e.g. each column)
   ! keeps its own four seeds, so streams can be advanced independently.

   use ccpp_kinds,      only: kind_phys
   use iso_fortran_env, only: int64

   implicit none
   private

   public :: kiss_rng_seed  ! one seed from a default integer
   public :: kiss_rng_next  ! advance a stream, return the next number

   integer(int64), parameter :: mask32 = 4294967295_int64   ! 2**32 - 1
   integer(int64), parameter :: two31  = 2147483648_int64   ! 2**31
   integer(int64), parameter :: two32  = 4294967296_int64   ! 2**32

contains

   function kiss_rng_seed(seed) result(s)
      ! The unsigned 32-bit value of a default-integer seed (the bits that
      ! kissvec's C code reads as uint32_t)
      !$acc routine seq
      integer, intent(in) :: seed
      integer(int64)      :: s

      s = iand(int(seed, int64), mask32)

   end function kiss_rng_seed

   function kiss_rng_next(s1, s2, s3, s4) result(ran)
      ! Advance the stream's seeds (from kiss_rng_seed or earlier calls) and
      ! return its next random number, in (0, 1)
      !$acc routine seq
      integer(int64), intent(inout) :: s1, s2, s3, s4
      real(kind_phys)               :: ran

      integer(int64) :: isum

      ! seed1 = 69069U * seed1 + 1327217885U
      s1 = iand(69069_int64*s1 + 1327217885_int64, mask32)
      ! seed2 ^= seed2 << 13; seed2 ^= seed2 >> 17; seed2 ^= seed2 << 5
      s2 = iand(ieor(s2, ishft(s2, 13)), mask32)
      s2 = ieor(s2, ishft(s2, -17))
      s2 = iand(ieor(s2, ishft(s2, 5)), mask32)
      ! seed3 = 18000U * (seed3 & 65535U) + (seed3 >> 16)
      s3 = 18000_int64*iand(s3, 65535_int64) + ishft(s3, -16)
      ! seed4 = 30903U * (seed4 & 65535U) + (seed4 >> 16)
      s4 = 30903_int64*iand(s4, 65535_int64) + ishft(s4, -16)
      ! (int32_t)(seed1 + seed2 + (seed3 << 16) + seed4) * 2.328306E-10 + 0.5
      isum = iand(s1 + s2 + iand(ishft(s3, 16), mask32) + s4, mask32)
      if (isum >= two31) then
         isum = isum - two32
      end if
      ran = real(isum, kind_phys) * 2.328306e-10_kind_phys + 0.5_kind_phys

   end function kiss_rng_next

end module atmos_phys_kiss_rng
