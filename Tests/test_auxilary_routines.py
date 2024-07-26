import unittest
import numpy as np
from auxilary_routines import CompoundClass, d
import unittest
import matplotlib.pyplot as plt
import numpy as np
from typing import List
from auxilary_routines import *

class TestAuxilaryRoutines(unittest.TestCase):

    def test_compound_class_organic(self):
        result = CompoundClass("C6H12O6")
        self.assertIn("org_", result)
        self.assertIn("HCO_", result)
        self.assertIn("HCO6", result)
        self.assertIn("nC=6", result)
        self.assertIn("nH=12", result)
        self.assertIn("nO=6", result)
        self.assertIn("nN=0", result)

    def test_compound_class_inorganic(self):
        result = CompoundClass("H2O")
        self.assertIn("inorg", result)
        self.assertIn("nC=0", result)
        self.assertIn("nH=2", result)
        self.assertIn("nO=1", result)
        self.assertIn("nN=0", result)

    def test_compound_class_13c(self):
        result = CompoundClass("13CC5H12")
        self.assertIn("org_", result)
        self.assertIn("HC_", result)
        self.assertIn("nC=6", result)
        self.assertIn("nH=12", result)
        self.assertIn("nO=0", result)
        self.assertIn("nN=0", result)

    def test_compound_class_nitrogen(self):
        result = CompoundClass("C2H5NO2")
        self.assertIn("org_", result)
        self.assertIn("orgN_", result)
        self.assertIn("orgN1", result)
        self.assertIn("NHCO2", result)
        self.assertIn("nC=2", result)
        self.assertIn("nH=5", result)
        self.assertIn("nO=2", result)
        self.assertIn("nN=1", result)

    def test_compound_class_sulfate(self):
        result = CompoundClass("C2H5SO4")
        self.assertIn("org_", result)
        self.assertIn("orgSO4", result)
        self.assertIn("nC=2", result)
        self.assertIn("nH=5", result)
        self.assertIn("nO=4", result)
        self.assertIn("nN=0", result)

    def test_d_numpy_array(self):
        arr = np.array([[1, 2, 3], [4, 5, 6]])
        result = d(arr)
        self.assertEqual(result, [2, 3])

    def test_d_list(self):
        lst = [1, 2, 3, 4]
        result = d(lst)
        self.assertEqual(result, [4])

    def test_d_empty(self):
        result = d(42)
        self.assertEqual(result, [])

    def test_DispCT(self):
        # Test that DispCT runs without errors
        try:
            DispCT(0)
        except Exception as e:
            self.fail(f"DispCT raised an exception: {e}")

    def test_drucksieb_comma(self):
        result = drucksieb(".", 1)
        expected = [
            [0, 0, 0, 0, 0, 0, 0, 0],
            [0, -2, -2, 0, 0, 0, 0, 0],
            [0, -2, -2, 0, 0, 0, 0, 0],
            [-2, -2, -2, -2, -2, 0, 0, 0]
        ]
        self.assertEqual(result, expected)

    def test_drucksieb_number(self):
        result = drucksieb("1", 1)
        expected = [
            [0, 0, 0, -2, -2, 0, 0, 0],
            [0, 0, -2, -2, -2, 0, 0, 0],
            [0, 0, 0, -2, -2, 0, 0, 0],
            [0, 0, 0, -2, -2, 0, 0, 0],
            [0, -2, -2, -2, -2, -2, 0, 0]
        ]
        self.assertEqual(result, expected)

    def test_drucksieb_multiple_digits(self):
        result = drucksieb("123", 1)
        self.assertEqual(len(result), 6)
        self.assertTrue(all(len(row) == 8 for row in result))

    def test_drucksieb_fontsize(self):
        result = drucksieb("1", 2)
        self.assertEqual(len(result), 10)
        self.assertEqual(len(result[0]), 16)

    def test_drucksieb_all_digits(self):
        digits = "0123456789"
        result = drucksieb(digits, 1)
        self.assertGreater(len(result), 6 * len(digits))

    def test_drucksieb_empty_string(self):
        result = drucksieb("", 1)
        self.assertEqual(result, [[0, 0, 0, 0, 0, 0, 0, 0]])

    def test_drucksieb_invalid_char(self):
        result = drucksieb("A", 1)
        self.assertEqual(result, [[0, 0, 0, 0, 0, 0, 0, 0]])

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.create_test_files()

    def tearDown(self):
        for file in os.listdir(self.temp_dir):
            os.remove(os.path.join(self.temp_dir, file))
        os.rmdir(self.temp_dir)

    def create_test_files(self):
        for i in range(3):
            filename = f'test_file_{i}.h5'
            filepath = os.path.join(self.temp_dir, filename)
            with h5py.File(filepath, 'w') as f:
                f.attrs['HDF5 File Creation Time'] = f'2023-05-{i+1:02d} 10:00:00'
                f.create_dataset('/FullSpectra/TofData', data=np.random.rand(100, 100))
                f.create_dataset('/SPECdata/Intensities', data=np.random.rand(50, 50))

    def test_FileList(self):
        result = FileList(None, self.temp_dir)
        self.assertIsInstance(result, dict)
        self.assertIn('path', result)
        self.assertIn('t09', result)
        self.assertIn('disp', result)
        self.assertEqual(len(result['path']), 3)
        self.assertEqual(len(result['t09']), 3)
        self.assertEqual(len(result['disp']), 3)

    def test_FileLst(self):
        # Create test files for FileLst
        for i in range(3):
            filename = f'test_ppb20_{i:02d}0101000000_UL.fdt'
            open(os.path.join(self.temp_dir, filename), 'w').close()

        result = FileLst(self.temp_dir)
        self.assertIsInstance(result, dict)
        self.assertIn('path', result)
        self.assertIn('disp', result)
        self.assertIn('startTimes', result)
        self.assertEqual(len(result['path']), 3)
        self.assertEqual(len(result['disp']), 3)
        self.assertEqual(len(result['startTimes']), 3)

    def test_starttimes(self):
        files = [
            '/path/to/ppb20230101000000_file.txt',
            '/path/to/ppb20230102000000_file.txt',
            '/path/to/ppb20230103000000_file.txt'
        ]
        result = starttimes(files, 'ppb20')
        self.assertEqual(len(result), 3)
        self.assertEqual(result, [20230101000000.0, 20230102000000.0, 20230103000000.0])

    def test_formula(self):
        lib = [
            [1, 1, 1],
            [1, 2, 3],  # C12
            [0, 1, 0],  # C13
            [2, 4, 6],  # H1
            [1, 2, 3],  # O16
            [0, 1, 0],  # O17
            [0, 0, 1],  # O18
            [1, 1, 2],  # N14
            [0, 1, 0],  # N15
            [0, 1, 0],  # Hpos
            [1, 0, 1],  # eneg
            [0, 1, 0],  # Cl35
            [0, 0, 1],  # Cl37
            [1, 0, 1],  # F
            [0, 1, 0],  # S32O4
            [0, 0, 1],  # S34O4
            [1, 0, 0],  # Si28
            [0, 1, 0],  # Si29
            [0, 0, 1],  # Si30
        ]
        result = formula(lib)
        self.assertEqual(len(result), 3)
        self.assertEqual(result[0], 'CH2OF+')
        self.assertEqual(result[1], '13CC2H4O217ONClS32O4SiH+')
        self.assertEqual(result[2], 'C3H6O318O2NF+')

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.create_test_files()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def create_test_files(self):
        # Create test HDF5 files
        for i in range(2):
            filename = f'test_file_{i}.h5'
            filepath = os.path.join(self.temp_dir, filename)
            with h5py.File(filepath, 'w') as f:
                f.attrs['HDF5 File Creation Time'] = f'2023-05-{i+1:02d} 10:00:00'
                f.create_dataset('/SPECdata/AverageSpec', data=np.random.rand(100))
                f.create_dataset('NbrWaveforms', data=1000)
                f.create_dataset('fullspectra/sampleinterval', data=0.1)
                f.create_dataset('RawData/HPTDC/PoissonCorrect', data=1.0)
                f.create_dataset('fullspectra/SINGLE_ION_SIGNAL', data=0.5)
                f.create_dataset('RawData/HPTDC/PoissonTdcDeadTime', data=0.01)
                f.attrs['COMPUTER_ID'] = 'test_computer'

    def test_get_eng_data_tof1000(self):
        event = {}
        file = os.path.join(self.temp_dir, 'test_file_0.h5')
        result = get_eng_data(event, file)
        
        self.assertIsInstance(result, dict)
        self.assertEqual(result['instrument'], 'TOF1000')
        self.assertIn('SumSpec', result)
        self.assertIn('cycles', result)
        self.assertIn('duration', result)
        self.assertIn('StartTime', result)

    def test_get_eng_data_tof8000(self):
        event = {}
        file = os.path.join(self.temp_dir, 'test_file_1.h5')
        result = get_eng_data(event, file)
        
        self.assertIsInstance(result, dict)
        self.assertEqual(result['instrument'], 'TOF8000')
        self.assertIn('SumSpec', result)
        self.assertIn('cycles', result)
        self.assertIn('duration', result)
        self.assertIn('StartTime', result)

    def test_get_eng_data_file_not_found(self):
        event = {}
        file = os.path.join(self.temp_dir, 'non_existent_file.h5')
        
        with self.assertRaises(FileNotFoundError):
            get_eng_data(event, file)

    def test_GetFileExport_load_from_files(self):
        event = {
            'Label_6': {'uVALUE': [-9999]},
            'Label_7': {'uVALUE': [-9999]},
            'base_1': {'uVALUE': [-9999]},
            'base_2': {'uVALUE': [-9999]},
            'base_3': {'uVALUE': [-9999]},
            'base_4': {'uVALUE': [-9999]},
            'base_5': {'uVALUE': [-9999]},
            'base_6': {'uVALUE': [-9999]},
            'base_10': {'uVALUE': [-9999]},
            'base_12': {'uVALUE': [-9999]},
            'Label_5': {'uVALUE': {'resolution': 1.0, 'peakshape': 'gaussian'}}
        }
        destfolder = self.temp_dir
        Name1 = 'test'

        # Create mock files
        os.makedirs(os.path.join(destfolder, 'Export/UL/IonList'), exist_ok=True)
        os.makedirs(os.path.join(destfolder, 'Export/FL/IonList'), exist_ok=True)
        os.makedirs(os.path.join(destfolder, 'Export/UL/cps/ocorr'), exist_ok=True)
        os.makedirs(os.path.join(destfolder, 'Export/UL/cps'), exist_ok=True)
        os.makedirs(os.path.join(destfolder, 'Export/FL/cps/ocorr'), exist_ok=True)
        os.makedirs(os.path.join(destfolder, 'Export/FL/cps'), exist_ok=True)
        os.makedirs(os.path.join(destfolder, 'Export/UL/ppb/ocorr'), exist_ok=True)
        os.makedirs(os.path.join(destfolder, 'Export/UL/ppb'), exist_ok=True)
        os.makedirs(os.path.join(destfolder, 'Export/FL/ppb/ocorr'), exist_ok=True)
        os.makedirs(os.path.join(destfolder, 'Export/FL/ppb'), exist_ok=True)

        for file_path in [
            f'Export/UL/IonList/MassIDs_{Name1}UL.csv',
            f'Export/FL/IonList/MassIDs_{Name1}FL.csv',
            f'Export/UL/cps/ocorr/corrcps{Name1}UL.fdt',
            f'Export/UL/cps/cps{Name1}UL.fdt',
            f'Export/FL/cps/ocorr/corrcps{Name1}FL.fdt',
            f'Export/FL/cps/cps{Name1}FL.fdt',
            f'Export/UL/ppb/ocorr/corrppb{Name1}UL.fdt',
            f'Export/UL/ppb/ppb{Name1}UL.fdt',
            f'Export/FL/ppb/ocorr/corrppb{Name1}FL.fdt',
            f'Export/FL/ppb/ppb{Name1}FL.fdt'
        ]:
            with open(os.path.join(destfolder, file_path), 'w') as f:
                f.write('1,2,3\n')

        result = GetFileExport(event, destfolder, Name1)

        self.assertIsInstance(result, dict)
        self.assertIn('massUL', result)
        self.assertIn('massFL', result)
        self.assertIn('corrcpsUL', result)
        self.assertIn('cpsUL', result)
        self.assertIn('corrcpsFL', result)
        self.assertIn('cpsFL', result)
        self.assertIn('corrppbUL', result)
        self.assertIn('ppbUL', result)
        self.assertIn('corrppbFL', result)
        self.assertIn('ppbFL', result)
        self.assertEqual(result['resolution'], 1.0)
        self.assertEqual(result['peakshape'], 'gaussian')

    def test_getPrimIons_h3o_mode(self):
        masses = np.array([19.018, 21.0221, 38.0326, 39.0327])
        cps = np.array([[100, 200, 300, 400],
                        [50, 100, 150, 200],
                        [10, 20, 30, 40],
                        [5, 10, 15, 20]])
        mode = 0

        result = getPrimIons(masses, cps, mode)

        self.assertIn('A', result)
        self.assertIn('B', result)
        self.assertEqual(len(result['A']), 4)
        self.assertEqual(len(result['B']), 4)

    def test_getPrimIons_no_mode(self):
        masses = np.array([30.994, 47.9966])
        cps = np.array([[100, 200],
                        [50, 100]])
        mode = 5

        result = getPrimIons(masses, cps, mode)

        self.assertIn('A', result)
        self.assertIn('B', result)
        self.assertEqual(len(result['A']), 2)
        self.assertEqual(len(result['B']), 2)

    def test_getPrimIons_o2_mode(self):
        masses = np.array([33.9935])
        cps = np.array([[100, 200]])
        mode = 9

        result = getPrimIons(masses, cps, mode)

        self.assertIn('A', result)
        self.assertIn('B', result)
        self.assertEqual(len(result['A']), 2)
        self.assertEqual(len(result['B']), 2)
        self.assertTrue(np.all(result['B'] == 0))

    @patch('auxilary_routines.open')
    def test_getpar_float_value(self, mock_open):
        mock_open.side_effect = [
            io.StringIO('C:/PTRwid/test_params.txt\n'),
            io.StringIO('_test_param=42.0\n')
        ]

        result = getpar('test_param')
        self.assertEqual(result, 42.0)

    @patch('auxilary_routines.open')
    def test_getpar_string_value(self, mock_open):
        mock_open.side_effect = [
            io.StringIO('C:/PTRwid/test_params.txt\n'),
            io.StringIO('_test_param=test_value\n')
        ]

        result = getpar('test_param')
        self.assertEqual(result, 'test_value')

    @patch('auxilary_routines.open')
    def test_getpar_not_found(self, mock_open):
        mock_open.side_effect = [
            io.StringIO('C:/PTRwid/test_params.txt\n'),
            io.StringIO('_other_param=42.0\n')
        ]

        result = getpar('test_param')
        self.assertIsNone(result)

    @patch('auxilary_routines.open')
    def test_getstr(self, mock_open):
        mock_open.side_effect = [
            io.StringIO('C:/PTRwid/test_params.txt\n'),
            io.StringIO('_test_param=test_value\n')
        ]

        result = getstr('test_param')
        self.assertEqual(result, 'test_value')

    @patch('auxilary_routines.open')
    def test_get_inst_id(self, mock_open):
        mock_open.return_value.__enter__.return_value.read.return_value = 'TEST_INSTRUMENT_ID'

        result = get_inst_id()
        self.assertEqual(result, 'TEST_INSTRUMENT_ID')


    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.transmission_file = os.path.join(self.temp_dir, 'TransmissionArchive.txt')
        with open(self.transmission_file, 'w') as f:
            f.write('_default=0,1,2,3,4,5\n')
            f.write('_campaign1=10,11,12,13,14,15\n')
        
        self.hdf5_file = os.path.join(self.temp_dir, 'test.h5')
        with h5py.File(self.hdf5_file, 'w') as f:
            f.create_dataset('/FullSpectra/TofData', data=np.random.rand(100, 100))
            f.create_dataset('/SPECdata/Intensities', data=np.random.rand(50, 50))

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_get_tr_existing_campaign(self):
        result = get_tr('campaign1', 5)
        self.assertEqual(result, [11.0, 12.0, 13.0, 14.0, 15.0])

    def test_get_tr_nonexistent_campaign(self):
        with self.assertLogs(level='WARNING') as cm:
            result = get_tr('nonexistent_campaign', 5)
        self.assertEqual(result, [1.0, 2.0, 3.0, 4.0, 5.0])
        self.assertIn('WARNING:root:UniqueCampaignName not found. Using default parameter', cm.output)

    def test_get_tr_time_warning(self):
        with self.assertLogs(level='WARNING') as cm:
            get_tr('campaign1', 100)
        self.assertIn('WARNING:root:nearest transmission more than 30 days away', cm.output)

    def test_add_tr(self):
        add_tr('new_campaign', [20, 21, 22, 23, 24, 25])
        with open(self.transmission_file, 'r') as f:
            content = f.read()
        self.assertIn('_new_campaign=20,21,22,23,24,25', content)

    def test_isodist(self):
        lib = [1, 0, 2, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0]
        result = isodist(lib)
        self.assertIn('masses', result)
        self.assertIn('dist', result)
        self.assertEqual(len(result['masses']), len(result['dist']))

    def test_LoadMassRange(self):
        result = LoadMassRange(self.hdf5_file, 10, 20, 1, 0, 1, 0.1)
        self.assertIsInstance(result, np.ndarray)

    def test_LoadMassRange_error(self):
        result = LoadMassRange('nonexistent.h5', 10, 20, 1, 0, 1, 0.1)
        self.assertEqual(result, -9999)

    def test_getSplit_small_dataset(self):
        result = getSplit(self.hdf5_file)
        self.assertEqual(result['split'], 1)

    def test_getSplit_large_dataset(self):
        with h5py.File(self.hdf5_file, 'w') as f:
            f.create_dataset('/FullSpectra/TofData', data=np.random.rand(1000, 1000))
        result = getSplit(self.hdf5_file)
        self.assertGreater(result['split'], 1)

    def test_LoadCycles(self):
        result = LoadCycles(self.hdf5_file, '/FullSpectra/TofData', [0, 0], [10, 10])
        self.assertEqual(len(result), 2)
        self.assertIsInstance(result[1], np.ndarray)

    def test_m2t_time_to_mass(self):
        result = m2t(5000, 1, 0, 0.5, 1e-10, time=True)
        self.assertLess(result, 5000)

    def test_m2t_mass_to_time(self):
        result = m2t(100, 1, 0, 0.5, 1e-10, mass=True)
        self.assertGreater(result, 100)

    def test_MakeCsv_2D_data(self):
        data = [['a', 'b', 'c'], ['1', '2', '3']]
        csv_file = os.path.join(self.temp_dir, 'test.csv')
        MakeCsv(csv_file, data)
        with open(csv_file, 'r') as f:
            content = f.read()
        self.assertEqual(content.strip(), 'a,b,c\n1,2,3')

    def test_MakeCsv_1D_data(self):
        data = ['a', 'b', 'c']
        csv_file = os.path.join(self.temp_dir, 'test.csv')
        MakeCsv(csv_file, data)
        with open(csv_file, 'r') as f:
            content = f.read()
        self.assertEqual(content.strip(), 'a,b,c')

    def test_read_csv(self):
        csv_file = os.path.join(self.temp_dir, 'test.csv')
        with open(csv_file, 'w') as f:
            f.write('1,2,3\n4,5,6')
        result = read_csv(csv_file)
        self.assertEqual(result, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

    def test_read_csv_semicolon(self):
        csv_file = os.path.join(self.temp_dir, 'test.csv')
        with open(csv_file, 'w') as f:
            f.write('1;2;3\n4;5;6')
        result = read_csv(csv_file, semicolon=True)
        self.assertEqual(result, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

    def test_read_csv_decimal_comma(self):
        csv_file = os.path.join(self.temp_dir, 'test.csv')
        with open(csv_file, 'w') as f:
            f.write('1,5;2,5;3,5\n4,5;5,5;6,5')
        result = read_csv(csv_file, semicolon=True, decimal_comma=True)
        self.assertEqual(result, [[1.5, 2.5, 3.5], [4.5, 5.5, 6.5]])

    def test_read_csv_str(self):
        csv_file = os.path.join(self.temp_dir, 'test.csv')
        with open(csv_file, 'w') as f:
            f.write('a,b,c\nd,e,f')
        result = read_csv_str(csv_file)
        self.assertEqual(result, [['a', 'd'], ['b', 'e'], ['c', 'f']])

    def test_read_csv_str_variable(self):
        csv_file = os.path.join(self.temp_dir, 'test.csv')
        with open(csv_file, 'w') as f:
            f.write('a,b,c\nd,e\nf,g,h,i')
        result = read_csv_str(csv_file, variable=True)
        self.assertEqual(result, [['a', 'b', 'c', ''], ['d', 'e', '', ''], ['f', 'g', 'h', 'i']])


    def test_read_index_file_valid(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.ind', delete=False) as temp_file:
            temp_file.write("_Textt0=10\n_Tex1=vector1\n_Dropp1=cond1\n_Text1=value1\n_Textt1=if1\n_Texttt1=else1\n_Dext15=ind_value1\n_Dextt15=value_2_add1")
            temp_file_path = temp_file.name

        result = read_index_file(temp_file_path)
        self.assertIsInstance(result, dict)
        self.assertEqual(result['interval'], '10')
        self.assertEqual(result['VectorInd'][0], 'vector1')
        self.assertEqual(result['CondInd'][0], 'cond1')
        self.assertEqual(result['ValueInd'][0], 'value1')
        self.assertEqual(result['ifInd'][0], 'if1')
        self.assertEqual(result['elseInd'][0], 'else1')
        self.assertEqual(result['IndValue'][0], 'ind_value1')
        self.assertEqual(result['Values2add'][0], 'value_2_add1')

        os.unlink(temp_file_path)

    def test_read_index_file_invalid(self):
        result = read_index_file('nonexistent_file.txt')
        self.assertEqual(result, {'VectorInd': -9999})

    def test_ReadFloat_valid_file(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
            temp_file.write("2\n10\n20\n30\n40\n50")
            temp_file_path = temp_file.name

        result = ReadFloat(temp_file_path)
        self.assertIsInstance(result, list)
        self.assertEqual(result, [30.0, 40.0, 50.0])

        os.unlink(temp_file_path)

    def test_ReadFloat_nonexistent_file(self):
        result = ReadFloat('nonexistent_file.txt')
        self.assertEqual(result, -9999)

    def test_ReArrange_1D_list(self):
        data = [1, 2, 3, 4, 5]
        result = ReArrange(data)
        self.assertEqual(result, [data])

    def test_starttimes_valid(self):
        files = ['ppb20230101000000_file.txt', 'ppb20230102000000_file.txt']
        result = starttimes(files, 'ppb')
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)
        self.assertAlmostEqual(result[0], 5114.0, places=1)
        self.assertAlmostEqual(result[1], 5115.0, places=1)

    def test_starttimes_invalid(self):
        files = ['invalid_file1.txt', 'invalid_file2.txt']
        result = starttimes(files, 'ppb')
        self.assertEqual(result, 0)

    def test_str2vec_valid(self):
        result = str2vec('1.0,2.0,3.0')
        self.assertEqual(result, [1.0, 2.0, 3.0])

    def test_str2vec_invalid(self):
        result = str2vec('invalid')
        self.assertEqual(result, -9999)

    def test_t09_format1(self):
        result = t09('01/01/2010, 12:00:00 PM')
        self.assertAlmostEqual(result, 365.5, places=1)

    def test_t09_format2(self):
        result = t09('01-Jan-10 12:00:00 PM')
        self.assertAlmostEqual(result, 365.5, places=1)

    def test_t09_format3(self):
        result = t09('01/01/2010 12h 00m 00s')
        self.assertAlmostEqual(result, 365.5, places=1)

    def test_t09_format4(self):
        result = t09('01.01.2010 12:00:00')
        self.assertAlmostEqual(result, 365.5, places=1)

    def test_t09_format5(self):
        result = t09('2010.01.01-12h00m00s')
        self.assertAlmostEqual(result, 365.5, places=1)

    def test_t09_format6(self):
        result = t09('2010-01-01 12:00:00', format6=True)
        self.assertAlmostEqual(result, 365.5, places=1)

    def test_t09str_float_input(self):
        result = t09str(365.5)
        self.assertEqual(result, '2010.01.01-12h00m00s')

    def test_t09str_list_input(self):
        result = t09str([2010, 1, 1, 12, 0, 0])
        self.assertEqual(result, '2010.01.01-12h00m00s')

if __name__ == '__main__':
    unittest.main()