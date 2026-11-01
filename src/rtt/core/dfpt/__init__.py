"""DFPT stage (team DFT): Quantum ESPRESSO electron-phonon workflow, a2F parsers, checks.

Main entry points::

    from rtt.core.dfpt import (EPhWorkflow, EPhSettings, PseudoTable, Launcher,
                           collect_qe_eph, read_two_column, parse_qe_a2f_dos,
                           parse_qe_lambda_file, parse_epw_a2f, epw_a2f_to_spectra)

Every Spectrum produced here has energies in kelvin and records in ``provenance`` the source
file, the discretisation rule, the dropped omega <= 0 weight and any clipping (ARCH section 5).
Each QE smearing gives one alternative Spectrum (ARCH S5); nothing is selected silently.
"""
from .attached import SpectrumCheckError, load_attached_spectra
from .band_coupling import (BandCoupling, EPWLambdaStates, band_coupling_from_epw_lambda,
                            band_coupling_from_spectra, epw_lambda_dos_average, epw_omega_max_meV,
                            parse_epw_lambda_aniso, parse_epw_lambda_fs, parse_epw_summed_lambda)
from .epw_inputs import collect_ph_for_epw, nscf_kpoints_card, write_epw_aniso_input
from .checks import ConsistencyCheck, Tolerances, check_epw_spectrum, check_qe_spectrum, spread_report
from .grids import EPhGrids, grid_from_spacing, select_grids
from .launcher import Launcher, RunRecord, run_command
from .parsers import (EPWA2F, IncompleteFileError, QEA2FDos, QELambdaEntry, atoms_to_spectrum, epw_a2f_to_spectra,
                      grid_to_atoms, parse_epw_a2f, parse_qe_a2f_dos, parse_qe_lambda_file,
                      mode_negative_fractions, qe_a2f_dos_to_spectrum, read_two_column, spectrum_lambda,
                      spectrum_omega_log)
from .pseudo import PseudoTable, read_upf_header
from .structure import CellData, parse_pw_final_structure, to_celldata
from .units import to_kelvin_factor
from .workflow import (SMEARING_ADMISSION, STEPS, EPhResult, EPhSettings, EPhWorkflow, WorkflowError,
                       admit_smearings,
                       build_fortify_shim, collect_qe_eph, min_frequency_cm1,
                       workflow_from_config)

__all__ = [
    "SpectrumCheckError", "load_attached_spectra", "IncompleteFileError",
    "BandCoupling", "EPWLambdaStates", "band_coupling_from_epw_lambda",
    "band_coupling_from_spectra", "epw_lambda_dos_average", "epw_omega_max_meV",
    "parse_epw_lambda_aniso", "parse_epw_lambda_fs", "parse_epw_summed_lambda",
    "collect_ph_for_epw", "nscf_kpoints_card", "write_epw_aniso_input",
    "ConsistencyCheck", "Tolerances", "check_epw_spectrum", "check_qe_spectrum", "spread_report",
    "EPhGrids", "grid_from_spacing", "select_grids", "Launcher", "RunRecord", "run_command",
    "EPWA2F", "QEA2FDos", "QELambdaEntry", "atoms_to_spectrum", "epw_a2f_to_spectra",
    "grid_to_atoms", "parse_epw_a2f", "parse_qe_a2f_dos", "parse_qe_lambda_file",
    "qe_a2f_dos_to_spectrum", "read_two_column", "spectrum_lambda", "spectrum_omega_log",
    "PseudoTable", "read_upf_header", "CellData", "parse_pw_final_structure", "to_celldata",
    "to_kelvin_factor", "STEPS", "EPhResult", "EPhSettings", "EPhWorkflow", "WorkflowError",
    "build_fortify_shim", "collect_qe_eph", "min_frequency_cm1", "workflow_from_config",
    "SMEARING_ADMISSION", "admit_smearings", "mode_negative_fractions",
]
