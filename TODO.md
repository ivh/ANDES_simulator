# ANDES E2E Simulator - TODO

## Not Yet Tested

- [ ] Other bands (Y, J, H, B, V, IZ, U) - expected to work like R-band
- [ ] PSF processing tool (`andes-sim psf-process`)
- [ ] HDF generation (`andes-sim generate-hdf`) - slow, low priority
- [ ] Error handling edge cases
- [ ] Batch processing modes

## Known Issues

- **PyEchelle multi-CPU bug**: Must use `max_cpu=1` (workaround applied)
- **Numba warnings**: Cosmetic only, can be ignored

## Optional Future Work

- Add pytest test suite
- Report PyEchelle multi-CPU bug upstream
- Performance benchmarking
