"""
Calculations provided by aiida_adamant.

Register calculations via the "aiida.calculations" entry point in setup.json.
"""
from __future__ import annotations

from typing import Optional

from aiida.common import datastructures, CalcInfo, CodeInfo
from aiida.common.folders import Folder
from aiida.engine import CalcJob
from aiida.orm import SinglefileData, StructureData, List
from aiida.plugins import DataFactory
from pymatgen.util.typing import PathLike

from aiida_adamant.utils.defaults import KgrnDefaults

KgrnInputData = DataFactory('adamant.kgrn_data')


class KgrnCalculation(CalcJob):
    """
    AiiDA calculation plugin wrapping the diff executable.

    Simple AiiDA plugin wrapper for 'diffing' two files.
    """
    @classmethod
    def define(cls, spec):
        """Define inputs and outputs of the calculation."""
        super().define(spec)

        spec.input_namespace('kgrn',
                             required=True,
                             help='Structure related kgrn input files')

        spec.input('kgrn.structure',
                   required=True,
                   valid_type=StructureData,
                   help='Structure that should be calculated')

        spec.input('kgrn.params',
                   required=True,
                   valid_type=KgrnInputData,
                   help='Kgrn calculation params')

        # Transfer Matrix
        spec.input('kgrn.transfer_matrix',
                   required=True,
                   valid_type=List,
                   help='List of structure constants')

        # Shape function
        spec.input('kgrn.shape_function',
                   required=True,
                   valid_type=SinglefileData,
                   help='Shape function of the structure')

        # Madelung constant
        spec.input('kgrn.madelung_matrix',
                   required=True,
                   valid_type=SinglefileData,
                   help='Madelung matrix of the structure')

        # Atomic config
        spec.input('kgrn.atom_cfg',
                   required=True,
                   valid_type=SinglefileData,
                   help='Atomic configuration for this structure')

        spec.input('metadata.options.input_filename',
                   valid_type=str,
                   default=KgrnDefaults.INPUT_FILENAME)

        spec.input('metadata.options.output_filename',
                   valid_type=str,
                   default=KgrnDefaults.OUTPUT_FILENAME)


        spec.exit_code(100,
                       'ERROR_MISSING_OUTPUT_FILES',
                       message='Calculation did not produce all expected '
                       'output files.')

        print("DEFINE WAS CALLED")

    def prepare_for_submission(self, folder: Folder) -> CalcInfo:
        """
        Create input files.

        :param folder: an `aiida.common.folders.Folder` where the plugin
        should temporarily place all files needed by
            the calculation.
        :return: `aiida.common.datastructures.CalcInfo` instance
        """


        with folder.open(self.options.input_filename, 'w',
                         encoding='utf8') as handle:

            handle.write(self.create_input_file_string())

        # Write shape function



        codeinfo = CodeInfo()
        codeinfo.code_uuid = self.inputs.code.uuid
        codeinfo.stdin_name = self.options.input_filename
        codeinfo.stdout_name = self.options.output_filename

        calcinfo = CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.retrieve_list = [self.options.output_filename]

        return calcinfo

    def create_input_file_string(self):

        lines = []

        lines += self._get_control_section()

        lines += self._get_scfp_section()

        lines.append(
            "**********************************************************************")
        lines.append(
            "Sort:  information for alloy:                                    "
            "    *")
        lines.append(
            "******************************SS-screeining*|***Magnetic "
            "structure ***")
        lines.append(
            "Symb  IQ  IT ITA NRM  CONC      a_scr b_scr |Teta    Phi    FXM  "
            "m(split)")

        for site in self.structure:

            component_index = 0

            for comp in site.species:
                site_index = site.properties['site_index']
                neq_site_index = site.properties['neq_site_index']

                alpha = site.species.screening_params[comp].alpha
                beta = site.species.screening_params[comp].beta
                model_param = site.species.magnetic_params[comp].magnetic_model
                init_mag_mom = site.species.magnetic_params[comp].init_mag_mom

                prefactor = [-1, 1] \
                    if site.species.magnetic_params[comp].is_paramagnetic \
                    else [1]

                for factor in prefactor:
                    component_index += 1
                    line = ""
                    line += f"{comp:<4s}{site_index:4d}{neq_site_index:4d}"
                    line += f"{component_index:4d}{1:4d}"
                    line += f"{site.species[comp] / len(prefactor):10.6f} "
                    line += f"{alpha:6.3f}"
                    line += f"{beta:6.3f}"
                    line += f"{0.0:8.4f}{0.0:8.4f}  {model_param:1s}"
                    line += f"{factor * init_mag_mom:9.4f}"
                    lines.append(line)

        lines.append(
            "**********************************************************************")
        lines.append("Spin-spiral wave vector:")
        lines.append("qx....={:^10.6f}qy....={:^10.6f}qz....={:^10.6f}".format(
            self.params['qx'], self.params['qy'], self.params['qz']
        ))

        lines.append(
            "**********************************************************************")
        lines.append(
            "Atom:  information for atomic calculation:                       "
            "    *")
        lines.append(
            "**********************************************************************")

        lines += self._get_atomic_section()

        return "\n".join([line.strip() for line in lines])


    def _get_control_section(self, rel_path: PathLike = None):

        lines = []

        lines.append(
            r"KGRN      HP..= 0   !                              16 Nov 00")

        lines.append("JOBNAM...={}".format(self.job_name))

        line = ""
        line += "MSGL.={:<2d}STRT.={:<2s}FUNC.={:<4s}".format(
            self.params['MSGL'],
            self.params['STRT'],
            self.params['FUNC'])

        line += "EXPAN={:<2d}FCD.={:<2s}GPM.={:<2s}FSM.={}".format(
            self.params['EXPAN'],
            self.params['FCD'],
            self.params['GPM'],
            self.params['FSM'])

        lines.append(line)

        lines.append("FOR001={}".format(
            self.config_files.transfer_matrix.get_string(rel_path)))
        lines.append("FOR002={}".format(
            self.config_files.madelung_matrix.get_string(rel_path)))
        lines.append("DIR003={}".format(self.output_dirs.ctrl_dir.get_string(
            rel_path, is_dir=True)))
        lines.append("DIR006={}".format(
            self.output_dirs.output_dir.get_string(rel_path, is_dir=True)))
        lines.append("DIR010={}".format(
            self.output_dirs.full_chd_dir.get_string(rel_path, is_dir=True)))
        lines.append("DIR011={}".format(""))
        lines.append("DIR021={}".format(""))
        lines.append("DIR022={}".format(
            self.config_files.shape_matrix.get_string(rel_path)))
        lines.append("FOR098={}".format(
            self.config_files.atom_config.get_string(rel_path)))

        lines.append(self.params['COMMENT'])

        lines.append("*" * 70)
        lines.append(
            "SCFP:  information for self-consistency procedure:               "
            "    *")
        lines.append("*" * 70)

        return lines

    def _get_energy_mesh_sector(self):

        lines = []

        line = f"ZMSH...={self.params['ZMSH']:>2s} "
        line += f"NZ1..={self.params['NZ1']:>3d} "
        line += f"NZ2..={self.params['NZ2']:>3d} "
        line += f"NZ3..={self.params['NZ3']:>3d} "
        line += f"NRES.={self.params['NRES']:>3d} "
        line += f"NZD..={self.params['NZD']:>3d} "

        lines.append(line)

        return lines

    def _get_scfp_section(self):

        lines = []

        lines.append(
            f"NITER.={self.params['NITER']:3d} NLIN.={self.params['NLIN']:3d} "
            f"NCPA.={self.params['NCPA']:3d} NPRN....={self.params['NPRN']}")

        line = f"FRC...={self.params['FRC']:>3s} "
        line += f"DOS..={self.params['DOS']:>3s} "
        line += f"OPS..={self.params['OPS']:>3s} "
        line += f"AFM..={self.params['AFM']:>3s} "
        line += f"CRT..={self.params['CRT']:>3s} "
        line += f"STMP..={self.params['STMP']:>2s} "

        lines.append(line)

        line = f"Lmaxh.={self.params['Lmaxh']:>3d} "
        line += f"Lmaxt={self.params['Lmaxt']:>3d} "
        line += f"NFI..={self.params['NFI']:>3d} "
        line += f"FIXG.={self.params['FIXG']:>3d} "
        line += f"SHF..={self.params['SHF']:>3d} "
        line += f"SOFC.={self.params['SOFC']:>3s} "

        lines.append(line)

        line = f"KMSH...={self.params['KMSH']:>2s} "
        line += f"IBZ..={self.params['IBZ']:>3d} "
        line += f"NKX..={self.params['NKX']:>3d} "
        line += f"NKY..={self.params['NKY']:>3d} "
        line += f"NKZ..={self.params['NKZ']:>3d} "
        line += f"FBZ..={self.params['FBZ']:>3s} "

        lines.append(line)

        lines += self._get_energy_mesh_sector()

        line = f"DEPTH..={self.params['DEPTH']:>7.3f} "
        line += f"IMAGZ.={self.params['IMAGZ']:>7.3f} "
        line += f"EPS...={self.params['EPS']:>7.3f} "
        line += f"ELIM..={self.params['ELIM']:>7.3f}"

        lines.append(line)

        line = f"AMIX...={self.params['AMIX']:>7.3f} "
        line += f"VMIX..={self.params['VMIX']:>7.3f} "
        line += f"EFMIX.={self.params['EFMIX']:>7.3f} "
        line += f"VMTZ..={self.params['VMTZ']:>7.3f}"

        lines.append(line)

        line = ""
        line += f"TOLE...={self.params['TOLE']:>6.1e} ".replace("e", "d")
        line += f"TOLEF.={self.params['TOLEF']:>6.1e} ".replace("e", "d")
        line += f"TOLCPA={self.params['TOLCPA']:>6.1e} ".replace("e", "d")
        line += f"TFERMI={self.params['TFERMI']:>7.1f} (K)"

        lines.append(line)

        line = f"SWS....={self.params['SWS']:>7.3f} "
        line += f"MMOM..={self.params['MMOM']:>7.3f}"

        lines.append(line)

        line = f"EFGS...={self.params['EFGS']:>7.3f} "
        line += f"HX....={self.params['HX']:>7.3f} "
        line += f"NX...={self.params['NX']:>3d} "
        line += f"NZ0..={self.params['NZ0']:>3d} "
        line += f"KPOLE={self.params['KPOLE']:>3d}"

        lines.append(line)

        return lines

    def _get_atomic_section(self):

        lines = []

        line = ""
        line += f"IEX...={self.params['IEX']:>3d} "
        line += f"NES..={self.params['NES']:>3d} "
        line += f"NITER={self.params['ANITER']:>3d} "
        line += f"IWAT.={self.params['IWAT']:>3d} "
        line += f"NPRNA={self.params['NPRNA']:>3d} "
        line += f"SOLV.={self.params['SOLV']:>3d} "

        lines.append(line)

        line = ""
        line += f"VMIXATM..={self.params['VMIXATM']:>10.6f} "
        line += f"RWAT....={self.params['RWAT']:>10.6f} "
        line += f"RMAX....={self.params['RMAX']:>10.6f} "
        lines.append(line)

        line = ""
        line += f"DX.......={self.params['DX']:>10.6f} "
        line += f"DR1.....={self.params['DR1']:>10.6f} "
        line += f"TEST....={self.params['TEST']:>10.2E} "
        lines.append(line)

        line = ""
        line += f"TESTE....={self.params['TESTE']:>10.2E} "
        line += f"TESTY...={self.params['TESTY']:>10.2E} "
        line += f"TESTV...={self.params['TESTV']:>10.2E} "
        lines.append(line)

        return lines




