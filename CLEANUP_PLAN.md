# ANDES E2E src Directory Cleanup Plan

**Status:** Phase 4 completed (Framework Validation - Flat Field)
**Date Started:** 2025-12-11
**Last Updated:** 2025-12-12

---

## Current State Analysis

The src directory contains:
- **andes_simulator/** - New consolidated framework (now tested and working)
- **15 legacy scripts** across 4 directories with ~60-70% code duplication
- **Data directories** (HDF/, SED/) with ZEMAX models and spectral data
- **Housekeeping issues** (log files, empty directories, naming inconsistencies)

## Cleanup Strategy

Since the new framework is **being validated**, we're taking a cautious approach:
1. Organize without breaking existing workflows ✅ DONE
2. Clean up obvious issues ✅ DONE
3. Prepare for framework validation ✅ DONE
4. Create clear migration path ✅ DONE
5. Validate framework ⏳ IN PROGRESS

---

## Phase 1: Immediate Housekeeping ✅ COMPLETED

### 1.1 Remove Temporary Files
**Files to remove:**
- `esorex.log` - temporary log file

**Action:** Delete and add to .gitignore

### 1.2 Handle Empty/Problematic Directories
**Directories:**
- `RIZ/` - empty, likely old output directory
- `RIZ_motecarlo/` - empty with typo ("motecarlo" → "montecarlo")

**Action:** Removed

### 1.3 Update .gitignore
**Add patterns:**
```
*.log
esorex.log
__pycache__/
*.pyc
.venv/
*.egg-info/
```

**Status:** ✅ Completed (Commit ccbd99c)

---

## Phase 2: Archive Legacy Scripts ✅ COMPLETED

### 2.1 Create Archive Structure
```
src/
├── legacy/
│   ├── README.md              # Migration guide pointing to andes_simulator
│   ├── FF_code/               # 6 flat field scripts
│   ├── FP_code/               # 2 Fabry-Perot scripts
│   ├── YJH/                   # 3 HDF generation scripts
│   ├── PSF/                   # 2 post-processing scripts
│   ├── pyechelle_test_ANDES.py
│   └── Rband_varyZemax.py
```

### 2.2 Legacy README Content
Document:
- Why scripts were archived
- Mapping to new framework commands
- How to run legacy scripts if needed
- Reference to andes_simulator/ORIGINAL_SCRIPT_MAPPING.md

### 2.3 Update Main Documentation
**Update CLAUDE.md** to indicate:
- Primary tool is now andes_simulator
- Legacy scripts available in legacy/ directory
- Framework testing status

**Status:** ✅ Completed (Commit 0cf1cd0)

---

## Phase 3: Documentation ✅ COMPLETED

### 3.1 Files Created/Updated
- `src/README.md` - Entry point for new users
- `src/legacy/README.md` - Archive explanation
- `src/CLAUDE.md` - Updated with new structure
- Framework documentation in andes_simulator/

**Status:** ✅ Completed (Commit 6661d90)

---

## Phase 4: Framework Validation ✅ COMPLETED (Flat Field)

### 4.1 Critical Fixes Applied

**Fix 1: PyEchelle API Compatibility** ✅
- Updated `Constant` → `ConstantPhotonFlux`
- Updated `CSV` → `CSVSource`
- Fixed all import statements
- Commit: 9142202

**Fix 2: Package Structure** ✅
- Moved `pyproject.toml` to `src/`
- Proper package installation with `uv sync`
- Commit: 9142202

**Fix 3: CLI Performance** ✅
- Lazy imports for fast startup
- `--help` now ~2s instead of several seconds
- Commit: 9cfd33c

**Fix 4: Flux Unit Handling** ✅ CRITICAL
- Problem: `ConstantPhotonFlux` expects ph/s/Å, not ph/s
- Solution: Added automatic flux unit conversion
- Supports: ph/s, ph/s/AA, ph/s/nm
- CLI defaults to ph/s/AA (no conversion warning)
- Commit: f7f91b3, 2ceb04c

**Fix 5: CPU Optimization** ✅
- Changed default from 1 core to 10 cores (performance cores only)
- Later reverted to max_cpu=1 due to PyEchelle multi-CPU bug (commit 9d61719)
- Single CPU provides acceptable performance (~18s for flat field)
- Commit: dc0ccc9 (optimization), 9d61719 (revert)

**Fix 6: Single Fiber Mode** ✅ CRITICAL
- Problem: Single fiber mode produced 0 photons
- Root cause: All fibers shared same dark source object
- Solution: Create individual `ConstantPhotonFlux` objects for each fiber
- PyEchelle requires separate object instances per fiber
- Commit: c48c4d2

### 4.2 Test Cases (Original Plan)

**Test 1: Flat Field - All Fiber Modes** ✅ COMPLETED
- Band tested: R-band (9232x9216, 66 fibers)
- Modes validated:
  - ✅ All fibers: 57M non-zero pixels, proper Poisson stats
  - ✅ Single fiber: 618k non-zero pixels (after Fix 6)
  - ✅ Even/odd: Two files created with correct fiber separation
  - ✅ Slit modes: first_slit (fibers 1-31), second_slit (fibers 35-66)
- **Performance:** ~18s with max_cpu=1
- **Parameters:** flux=100 ph/s/AA, exposure=10-50s
- **See:** VALIDATION_REPORT.md for complete details

**Test 2: Fabry-Perot Calibration** ✅ COMPLETED
- ✅ Fixed flux units: "ph/s" instead of "ph/s/AA"
- ✅ Fixed scaling: pre-process CSV data before CSVSource creation
- ✅ Fixed object sharing: individual sources per fiber
- ✅ R-band validated: 66,877 photons with flux=100, proper statistics
- ✅ Performance: 2.8s for single fiber
- ✅ Added user-friendly `--flux` parameter

**Test 3: Post-Processing** ⏸️ DEFERRED
- Requires Fabry-Perot outputs
- Will test after FP issues resolved

**Test 4: HDF Generation** ⏸️ DEFERRED (Non-critical, very slow)
- Test on one band only if needed
- Verify ZEMAX interface works correctly

### 4.3 Validation Checklist (Flat Field Mode)

- [x] Framework runs without errors
- [x] Output files created with correct dimensions (9232x9216 for R-band)
- [x] Non-zero data produced (proper Poisson statistics)
- [x] Flux units handled correctly (ph/s/AA)
- [x] Performance acceptable (~18s for R-band with max_cpu=1)
- [x] Documentation clear and updated
- [x] R-band fully validated (other bands expected to work similarly)
- [x] All fiber illumination modes tested (all, single, even/odd, slits)
- [ ] Fabry-Perot mode - deferred (requires investigation)
- [ ] Post-processing - deferred (needs FP outputs)
- [ ] Error handling - not tested
- [ ] Batch processing - not tested

### 4.4 Known Issues

**Issue 1: PyEchelle Multi-CPU Bug** ⚠️ CRITICAL
- PyEchelle 0.4.0 has bug in multi-CPU result aggregation
- Error: `TypeError: 'int' object is not subscriptable` in `_simulate_multi_cpu`
- Occurs at line 364: `ccd_results = [r[0] for r in results]`
- **Workaround:** Default max_cpu=1 (commit 9d61719)
- **Impact:** Single CPU provides acceptable performance (~18s)
- **Status:** Workaround applied, may report to PyEchelle maintainers

**Issue 2: Legacy Script Incompatibility**
- Legacy scripts use PyEchelle <0.4.0 API
- Cannot run with current PyEchelle 0.4.0
- **Resolution:** Accept functional validation instead of numerical comparison
- **Status:** Documented, accepted

**Issue 3: Numba Threading Warnings**
- Workqueue threading layer warnings
- Semaphore cleanup warnings
- **Impact:** Cosmetic only, no functional impact
- **Resolution:** TBB not available for macOS ARM; warnings can be ignored
- **Status:** Documented, accepted

**Issue 4: Fabry-Perot CSV Source** ⏸️ INVESTIGATING
- FP mode produces 0 photons despite correct CSV data
- CSV file exists and covers wavelength range
- Scaling factors tested extensively
- Possible PyEchelle 0.4.0 CSV source API issue
- **Status:** Requires separate investigation, not blocking flat-field use

---

## Phase 5: Final Structure (After Validation) ⏸️ PENDING

### 5.1 Target Directory Structure
```
src/
├── andes_simulator/          # Main framework (production-ready)
│   ├── core/                 # Simulation engine
│   ├── sources/              # Source models
│   ├── models/               # HDF builders
│   ├── postprocess/          # Post-processing
│   ├── cli/                  # Command interface
│   ├── configs/              # Configuration templates
│   └── README.md             # Usage documentation
├── HDF/                      # ZEMAX optical models
├── SED/                      # Spectral energy distributions
├── legacy/                   # Archived scripts (reference only)
├── CLAUDE.md                 # Project context
├── CLEANUP_PLAN.md           # This file
├── VALIDATION_REPORT.md      # Validation findings
├── .gitignore                # Ignore patterns
└── README.md                 # High-level overview
```

### 5.2 Documentation Updates
**Create/update:**
- [x] `src/README.md` - Entry point for new users
- [x] `src/legacy/README.md` - Archive explanation
- [x] Framework status in CLAUDE.md
- [x] `VALIDATION_REPORT.md` - Detailed validation results
- [ ] Quick start guide for common tasks (can be added later)

---

## Phase 6: Optional Optimizations (Future) ⏸️ DEFERRED

### 6.1 Consolidate Duplicate Code (If Keeping Legacy)
If legacy scripts must remain active:
- Extract common code to `legacy/common.py`
- Create unified scripts with mode parameters
- Reduce 15 scripts → 5-6 scripts

**Status:** Not needed - framework is working

### 6.2 Enhance Framework
- Add validation mode that compares outputs to legacy
- Performance profiling
- Extended test coverage
- CI/CD integration

### 6.3 Data Organization
Keep HDF/ and SED/ as-is per user preference.

---

## Critical Files Modified

**Configuration:**
- `src/.gitignore` - Created
- `src/pyproject.toml` - Moved from andes_simulator/
- `src/uv.lock` - Moved from andes_simulator/

**Framework Code:**
- `andes_simulator/core/simulator.py` - Added flux unit conversion
- `andes_simulator/core/config.py` - Updated max_cpu default to 10
- `andes_simulator/cli/main.py` - Lazy imports, flux_unit fix
- `andes_simulator/sources/*.py` - Updated PyEchelle API imports
- `andes_simulator/configs/examples/*.yaml` - Updated flux units

**Documentation:**
- `src/README.md` - Created
- `src/CLAUDE.md` - Updated
- `src/legacy/README.md` - Created
- `src/VALIDATION_REPORT.md` - Created
- `src/CLEANUP_PLAN.md` - This file

**Files Moved:**
- `FF_code/` → `legacy/FF_code/`
- `FP_code/` → `legacy/FP_code/`
- `YJH/` → `legacy/YJH/`
- `PSF/` → `legacy/PSF/`
- `pyechelle_test_ANDES.py` → `legacy/`
- `Rband_varyZemax.py` → `legacy/`

**Files Deleted:**
- `esorex.log`
- `RIZ/` (empty directory)
- `RIZ_motecarlo/` (empty directory)

---

## Recommended Execution Order

1. ✅ **Phase 1** - Safe housekeeping (no risk)
2. ✅ **Phase 2** - Archive legacy scripts (reversible via git)
3. ✅ **Phase 3** - Create documentation
4. ⏳ **Phase 4** - Validate framework (IN PROGRESS)
5. ⏸️ **Phase 5** - Final polish (after validation complete)
6. ⏸️ **Phase 6** - Optional enhancements (future)

---

## Risk Assessment

**Low Risk:**
- Removing log files ✅ Done
- Updating .gitignore ✅ Done
- Creating documentation ✅ Done

**Medium Risk:**
- Moving legacy scripts to archive ✅ Done (reversible via git)
- Removing empty directories ✅ Done

**High Risk:**
- None (framework validation happened before any breaking changes)

---

## Success Criteria

- [x] No existing workflows are broken
- [x] Directory structure is clear and logical
- [x] New framework runs successfully
- [x] Legacy scripts remain accessible in archive
- [x] Documentation guides users to new framework
- [x] Repository is clean and professional
- [x] Critical bugs fixed (flux units, API compatibility)
- [ ] Comprehensive validation complete (partial)

---

## Git Commit History

```
2ceb04c Update example config to use correct flux units
dc0ccc9 Optimize CPU usage: use 10 performance cores by default
d369f9e Update validation report with flux unit fix results
f7f91b3 Fix flux unit handling for ConstantPhotonFlux sources
f5e8e01 Phase 4: Document initial validation findings
9cfd33c Optimize CLI startup with lazy imports
9142202 Phase 4 prep: Fix pyechelle API compatibility
6661d90 Phase 3: Add comprehensive documentation
0cf1cd0 Phase 2: Archive legacy scripts to legacy/ directory
ccbd99c Phase 1: Clean up temporary files and add .gitignore
874376f Pre-cleanup snapshot: current state of src directory
```

---

## Next Steps

1. **Complete Phase 4 validation** (in progress)
   - Test Fabry-Perot simulations
   - Test single fiber modes
   - Test post-processing tools
   - Document any additional issues

2. **Update VALIDATION_REPORT.md** with complete findings

3. **Phase 5: Final polish** (when validation complete)
   - Any remaining documentation updates
   - Performance benchmarking
   - User migration guide if needed

4. **Push to remote** when ready to share

---

*Plan created: 2025-12-11*
*Currently in Phase 4 (Framework Validation)*
*Framework Status: ✅ Core functionality working, validation in progress*
