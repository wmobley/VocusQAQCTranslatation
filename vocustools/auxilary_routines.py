# Standard library imports
import os
import re
from datetime import datetime, timedelta
from typing import List, Union, Dict, Tuple, Any, Optional
import time
import glob
# Third-party imports
import h5py
import matplotlib.pyplot as plt
import numpy as np
# Configuration
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['font.size'] = 12


def CompoundClass(formula: str) -> List[str]:
    # Lines 11857-11866: Initialize variables
    C, S, N, O, H = 0, 0, 0, 0, 0
    
    # Lines 11868-11876: Check for 13CC
    if '13CC' in formula:
        numb1 = formula[formula.index('13CC') + 4]
        numb2 = formula[formula.index('13CC') + 5]
        if numb1.isdigit() and numb2.isdigit():
            C = int(numb1) * 10 + int(numb2) + 1
        elif numb1.isdigit():
            C = int(numb1) + 1
        else:
            C = 2

    # Lines 11878-11885: Check for C (not 13CC)
    elif 'C' in formula and '13CC' not in formula:
        numb1 = formula[formula.index('C') + 1]
        numb2 = formula[formula.index('C') + 2]
        if numb1.isdigit() and numb2.isdigit():
            C = int(numb1) * 10 + int(numb2)
        elif numb1.isdigit():
            C = int(numb1)
        else:
            C = 1

    # Line 11887: If no C, set C to 0
    if 'C' not in formula:
        C = 0

    # Line 11889: Check for SO4
    S = 1 if 'SO4' in formula else 0

    # Lines 11891-11895: Check for N
    if 'N' in formula:
        numb = formula[formula.index('N') + 1]
        N = int(numb) if numb.isdigit() else 1
    else:
        N = 0

    # Lines 11897-11901: Check for O
    if 'O' in formula:
        numb = formula[formula.index('O') + 1]
        O = int(numb) if numb.isdigit() else 1
    else:
        O = 0

    # Lines 11903-11911: Check for H
    if 'H' in formula:
        numb1 = formula[formula.index('H') + 1]
        numb2 = formula[formula.index('H') + 2]
        if numb1.isdigit() and numb2.isdigit():
            H = int(numb1) * 10 + int(numb2)
        elif numb1.isdigit():
            H = int(numb1)
        elif numb1 != '+':
            H = 1
        else:
            H = 0

    # Lines 11913-11942: Determine compound type
    type_list = []
    if C == 0:
        type_list.append('inorg')
    if C >= 1:
        type_list.append('org_')
    if C >= 1 and O == 0 and N == 0 and S == 0:
        type_list.append('HC_')
    if C >= 1 and S >= 1:
        type_list.append('orgSO4')
    if C >= 1 and S == 0 and N >= 1:
        type_list.append('orgN_')
    if C >= 1 and S == 0 and N == 1:
        type_list.append('orgN1')
    if C >= 1 and S == 0 and N == 2:
        type_list.append('orgN2')
    if C >= 1 and S == 0 and O == 0 and N >= 1:
        type_list.append('NHCO0')
    for i in range(1, 9):
        if C >= 1 and S == 0 and O == i and N >= 1:
            type_list.append(f'NHCO{i}')
    if C >= 1 and S == 0 and N == 0:
        type_list.append('orgNoN')
    if C >= 1 and S == 0 and O >= 1 and N == 0:
        type_list.append('HCO_')
    for i in range(1, 9):
        if C >= 1 and S == 0 and O == i and N == 0:
            type_list.append(f'HCO{i}')

    # Lines 11944-11949: Add additional information
    unsat = f"unsat={((2.0 + 2*float(C) + float(N) - float(H))/2):.1f}{'N' if N > 0 else 'C'}"
    type_list.append(unsat)

    N_ox_st = -1  # This should be replaced with the actual value from getpar('N_ox_st')
    OSC = f"OSC={((2*float(O)/float(C) - N_ox_st*float(N)/float(C) - float(H)/float(C))):.5f}"
    type_list.append(OSC)

    type_list.extend([
        f"nC={C:2d}",
        f"nH={H:2d}",
        f"nO={O:2d}",
        f"nN={N:2d}"
    ])

    return type_list[1:]  # Return all elements except the first (empty) one

def d(array: Union[List, np.ndarray]) -> List[int]:
    """
    # Line 11947-11949: Function to get dimensions of an array
    Equivalent to IDL's size(array, /dimensions)
    """
    if isinstance(array, np.ndarray):
        return list(array.shape)
    elif isinstance(array, list):
        return [len(array)]
    else:
        return []

def DispCT(index: int) -> None:
    """
    # Lines 11953-11962: Display color table
    """
    
    # Line 11954: Set up multi-plot layout
    fig, axs = plt.subplots(6, 1, figsize=(10, 20))
    
    # Line 11955: Set color mode (not needed in matplotlib)
    
    # Line 11956: Load color table (using matplotlib's color maps)
    cmap = plt.get_cmap('viridis')  # Example colormap, adjust as needed
    
    # Lines 11957-11961: Create color display
    for k in range(6):
        ax = axs[k]
        for i in range(-25, 26):
            color = cmap((i + 25 + k * 50) / 300)
            ax.axvline(i + k * 50, color=color, linewidth=15)
        ax.set_xlim(-25 + k * 50, 25 + k * 50)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
    
    plt.tight_layout()
    plt.show()

def drucksieb(string: str, fontsize: int) -> List[List[int]]:
    """
    Convert a string of numbers into a matrix representation for image writing.
    
    Args:
    string (str): The input string to convert
    fontsize (int): The size of the font to use
    
    Returns:
    List[List[int]]: A 2D matrix representation of the string
    """
    # Lines 11965-11975: Define digit patterns
    komma = [1, 2, 4, 5]
    eins = [3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 26, 27, 28, 29, 33, 34, 39]
    zwei = [1, 2, 3, 4, 5, 7, 8, 14, 15, 21, 22, 28, 34, 35, 37, 38, 40, 41, 44, 45, 46]
    drei = [2, 3, 4, 7, 8, 10, 11, 16, 17, 21, 22, 27, 28, 34, 35, 37, 38, 40, 41, 44, 45, 46]
    vier = [4, 10, 13, 14, 15, 16, 17, 19, 22, 26, 28, 32, 34, 39, 40, 45, 46]
    fuenf = [2, 3, 4, 7, 8, 10, 11, 16, 17, 22, 23, 25, 26, 27, 28, 31, 32, 37, 38, 43, 44, 45, 46, 47]
    sechs = [3, 4, 8, 10, 11, 13, 14, 16, 17, 19, 20, 22, 23, 25, 26, 27, 28, 31, 32, 38, 40, 45, 46]
    sieben = [2, 3, 8, 9, 15, 21, 27, 28, 33, 34, 40, 41, 43, 44, 45, 46, 47]
    acht = [2, 3, 4, 7, 8, 10, 11, 13, 14, 16, 17, 20, 21, 22, 26, 27, 28, 31, 32, 34, 35, 37, 38, 40, 41, 44, 45, 46]
    neun = [2, 3, 8, 10, 16, 17, 20, 21, 22, 23, 25, 26, 28, 29, 31, 32, 34, 35, 37, 38, 40, 44, 45]
    null = [3, 8, 10, 13, 14, 16, 17, 19, 20, 22, 23, 25, 26, 28, 29, 31, 32, 34, 35, 38, 40, 45]

    # Line 11977: Get length of input string
    length = len(string)

    # Line 11978: Initialize siebmat
    siebmat: List[List[int]] = [[0] * 8]

    # Lines 11979-12003: Process each character in the string
    for i in range(length):
        substri = string[i]
        submat: List[List[int]] = [[0] * 8 for _ in range(6)]

        if substri == '.':
            submat = [[0] * 8 for _ in range(3)]
            for idx in komma:
                submat[idx // 3][idx % 3] = -2
        elif substri == '1':
            submat = [[0] * 8 for _ in range(5)]
            for idx in eins:
                submat[idx // 5][idx % 5] = -2
        elif substri == '2':
            for idx in zwei:
                submat[idx // 6][idx % 6] = -2
        elif substri == '3':
            for idx in drei:
                submat[idx // 6][idx % 6] = -2
        elif substri == '4':
            for idx in vier:
                submat[idx // 6][idx % 6] = -2
        elif substri == '5':
            for idx in fuenf:
                submat[idx // 6][idx % 6] = -2
        elif substri == '6':
            for idx in sechs:
                submat[idx // 6][idx % 6] = -2
        elif substri == '7':
            for idx in sieben:
                submat[idx // 6][idx % 6] = -2
        elif substri == '8':
            for idx in acht:
                submat[idx // 6][idx % 6] = -2
        elif substri == '9':
            for idx in neun:
                submat[idx // 6][idx % 6] = -2
        elif substri == '0':
            for idx in null:
                submat[idx // 6][idx % 6] = -2

        siebmat.extend(submat)

    # Lines 12004-12005: Get dimensions of siebmat
    dim = [len(siebmat), len(siebmat[0])]

    # Line 12006: Expand siebmat
    siebmatrix = [[siebmat[i // fontsize][j // fontsize] for j in range(dim[1] * fontsize)] for i in range(dim[0] * fontsize)]

    # Line 12007: Set first 5 elements to -2 if all elements are greater than -0.5
    if all(all(val > -0.5 for val in row) for row in siebmatrix):
        for i in range(5):
            siebmatrix[0][i] = -2

    # Line 12008: Set all elements less than -0.5 to -9999
    siebmatrix = [[val if val >= -0.5 else -9999 for val in row] for row in siebmatrix]

    return siebmatrix

def FileList(event, Path: str) -> Dict[str, List]:
    # Line 12010-12011: Function definition and compilation option
    # IDL's compile_opt idl2 is not needed in Python

    # Line 12013-12017: Initialize variables
    Err: List[str] = ['']
    locind: int = 0
    addind: int = 0
    location: str = '/FullSpectra/TofData'
    oldfile: str = '-9999'
    mist: float = time.time()

    # Line 12018-12019: Search for files
    Recur_Pattern: str = '*.h5'
    Files: List[str] = [os.path.join(root, name)
                        for root, dirs, files in os.walk(Path)
                        for name in files
                        if name.endswith('.h5')]

    # Line 12021-12022: Get number of files
    n: int = len(Files)

    # Line 12036-12037: Initialize error tracking
    errFiles: List[int] = [-1]
    jj: int = 0

    # Line 12039-12040: Set initial location
    location: str = '/FullSpectra/TofData'

    # Line 12041-12049: Error handling setup (Python equivalent of IDL's CATCH)
    def handle_error(error_status: int) -> None:
        nonlocal location, errFiles, jj, Err
        if location == '/SPECdata/Intensities':
            location = '/FullSpectra/TofData'
            errFiles.append(jj)
        else:
            location = '/SPECdata/Intensities'
            errFiles.append(jj - 1)
        print('')

    # Line 12051-12052: Update error log
    if max(errFiles) == jj:
        Err.extend(['', f'Discard: {Files[max(errFiles)]}'])

    # Line 12053-12076: Main file processing loop
    ind: List[int] = []
    for jj in range(max(errFiles) + 1, n):
        if min(abs(np.array(errFiles) - jj)) > 0:
            try:
                with h5py.File(Files[jj], 'r') as file_id:
                    dataset_id1 = file_id[location]
                    Dims = dataset_id1.shape
                    Start: int = 10000
                    Width: int = 2
                    if len(Dims) >= 2:
                        ddd = [0] * (len(Dims) - 1)
                        eee = list(Dims[1:])
                        start2 = [Start] + ddd
                        count = [Width] + [d + 1 for d in ddd]
                    else:
                        start2 = [Start]
                        count = [Width]
                    ptrData = dataset_id1[tuple(slice(s, s+c) for s, c in zip(start2, count))]
                ind.append(jj)
                # Update info widget (Python equivalent needed)
                # WIDGET_CONTROL(WIDGET_INFO(event.TOP, FIND_BY_UNAME='Text_info'), SET_VALUE=[info, Err])
            except Exception as e:
                handle_error(1)

    # Line 12080-12081: Update Files list
    Files = [Files[i] for i in ind] if ind else []
    n = len(Files)

    # Line 12083-12086: Initialize time-related variables
    CreationTime: List[str] = ['nix']
    StartTimes: List[float] = [0]
    StartTimesStr: List[str] = ['nix']

    # Line 12087-12108: Process file creation times
    if Files:
        for i, file in enumerate(Files):
            with h5py.File(file, 'r') as fid:
                mistt = 'nix'
                for attr_name, attr_value in fid.attrs.items():
                    if 'HDF5 File Creation Time' in attr_name:
                        mistt = attr_value
                        break
                if mistt == 'nix':
                    for attr_name, attr_value in fid.attrs.items():
                        if 'FileCreatedTimeSTR' in attr_name or 'FileCreatedTimeString' in attr_name:
                            mistt = attr_value
                            break
            CreationTime.append(mistt)
            tt09 = t09(mistt)  # Assuming t09 function is defined elsewhere
            StartTimes.append(tt09)
            StartTimesStr.append(t09str(tt09))  # Assuming t09str function is defined elsewhere
            # Update info widget (Python equivalent needed)
            # WIDGET_CONTROL(WIDGET_INFO(event.TOP, FIND_BY_UNAME='Text_info'), SET_VALUE=[info, info2, Err])

    # Line 12110: Remove first element of CreationTime
    CreationTime = CreationTime[1:]

    # Line 12114-12131: Sort files based on start times
    sorted_indices = sorted(range(len(StartTimes)), key=lambda k: StartTimes[k])
    Files = [Files[i] for i in sorted_indices]
    StartTimesStr = [StartTimesStr[i] for i in sorted_indices]
    StartTimes = [StartTimes[i] for i in sorted_indices]

    # Line 12135-12137: Create and return result structure
    s1: Dict[str, List] = {
        'path': Files,
        't09': StartTimes,
        'disp': StartTimesStr
    }

    return s1

def FileLst(Path: str) -> dict:
    """
    Retrieve file list and process it.
    
    Args:
    Path (str): The directory path to search for files.
    
    Returns:
    dict: A dictionary containing file paths, display names, and start times.
    """
    # Line 12146: Set recursive pattern
    Recur_Pattern = '*UL.fdt'
    
    # Line 12147: Search for files
    Files = glob.glob(os.path.join(Path, Recur_Pattern))
    
    # Lines 12148-12149: Filter files
    Files = [f for f in Files if 'ppb20' in f and 'corrppb20' not in f]
    
    # Lines 12151-12153: Sort files by start times
    StartTimes = starttimes(Files, 'ppb20')
    Files = [f for _, f in sorted(zip(StartTimes, Files))]
    StartTimes.sort()
    
    # Lines 12154-12157: Create display names
    mist3: List[str] = []
    for i, File in enumerate(Files):
        ppb_index = File.index('ppb20')
        mist3.append(f"{i:4d} - {File[ppb_index+3:ppb_index+23]}")
    
    # Lines 12159-12161: Final filtering
    filter_indices = [i for i, f in enumerate(Files) if 'UL.fdt' in f]
    Files = [Files[i] for i in filter_indices]
    mist3 = [mist3[i] for i in filter_indices]
    
    # Line 12162: Create return structure
    s1 = {
        'path': Files,
        'disp': mist3,
        'startTimes': StartTimes
    }
    
    # Line 12163: Return result
    return s1

def starttimes(files: List[str], prefix: str) -> List[float]:
    """
    Extract start times from file names.
    
    Args:
    files (List[str]): List of file paths.
    prefix (str): Prefix to look for in file names.
    
    Returns:
    List[float]: List of start times.
    """
    times: List[float] = []
    for file in files:
        try:
            # Extract time string after the prefix
            time_str = file.split(prefix)[1][:14]
            # Convert to float (assuming format is YYYYMMDDhhmmss)
            time_float = float(time_str)
            times.append(time_float)
        except (IndexError, ValueError):
            # If extraction fails, append a default value
            times.append(0.0)
    return times

def formula(lib: List[List[float]]) -> List[str]:
    # Lines 12183-12184: Determine the length of the input
    length = len(lib[0]) if len(lib) > 1 else 1
    
    # Lines 12185-12200: Initialize string arrays for each element
    C12 = [''] * length
    C13 = [''] * length
    H1 = [''] * length
    O16 = [''] * length
    O17 = [''] * length
    O18 = [''] * length
    N14 = [''] * length
    N15 = [''] * length
    Hpos = [''] * length
    eneg = [''] * length
    Cl35 = [''] * length
    Cl37 = [''] * length
    F = [''] * length
    S32O4 = [''] * length
    S34O4 = [''] * length
    Si = [''] * length

    # Line 12202: Get the number of rows in lib
    len_lib = len(lib)

    # Lines 12204-12206: Process C12
    FIL = [i for i, x in enumerate(lib[1]) if x > 0]
    if FIL:
        for i in FIL:
            C12[i] = 'C'
        for i in [i for i in FIL if lib[1][i] > 9.5]:
            C12[i] += f"{int(lib[1][i]):02d}"
        for i in [i for i in FIL if 1.5 < lib[1][i] < 9.5]:
            C12[i] += f"{int(lib[1][i]):1d}"

    # Lines 12208-12209: Process C13
    FIL = [i for i, x in enumerate(lib[2]) if x > 0]
    if FIL:
        for i in FIL:
            C13[i] = '13C'

    # Lines 12211-12214: Process H1
    FIL = [i for i, x in enumerate(lib[3]) if x > 0]
    if FIL:
        for i in FIL:
            H1[i] = 'H'
        for i in [i for i in FIL if lib[3][i] > 9.5]:
            H1[i] += f"{int(lib[3][i]):02d}"
        for i in [i for i in FIL if 1.5 < lib[3][i] < 9.5]:
            H1[i] += f"{int(lib[3][i]):1d}"

    # Lines 12216-12219: Process O16
    FIL = [i for i, x in enumerate(lib[4]) if x > 0]
    if FIL:
        for i in FIL:
            O16[i] = 'O'
        for i in [i for i in FIL if lib[4][i] > 9.5]:
            O16[i] += f"{int(lib[4][i]):02d}"
        for i in [i for i in FIL if 1.5 < lib[4][i] < 9.5]:
            O16[i] += f"{int(lib[4][i]):1d}"

    # Lines 12221-12222: Process O17 and O18
    for i, x in enumerate(lib[5]):
        if x > 0:
            O17[i] = '17O'
    for i, x in enumerate(lib[6]):
        if x > 0:
            O18[i] = '18O'

    # Lines 12225-12228: Process N14
    FIL = [i for i, x in enumerate(lib[7]) if x > 0]
    if FIL:
        for i in FIL:
            N14[i] = 'N'
        for i in [i for i in FIL if lib[7][i] > 9.5]:
            N14[i] += f"{int(lib[7][i]):02d}"
        for i in [i for i in FIL if 1.5 < lib[7][i] < 9.5]:
            N14[i] += f"{int(lib[7][i]):1d}"

    # Lines 12230-12231: Process N15
    for i, x in enumerate(lib[8]):
        if x > 0:
            N15[i] = '15N'

    # Lines 12233-12253: Process additional elements if len_lib > 11
    if len_lib > 11:
        # Process Cl35, Cl37, F, S32O4, S34O4, Si
        for elem, idx in [(Cl35, 11), (Cl37, 12), (F, 13), (S32O4, 14), (S34O4, 15)]:
            FIL = [i for i, x in enumerate(lib[idx]) if x > 0]
            if FIL:
                for i in FIL:
                    elem[i] = elem.__name__
                for i in [i for i in FIL if lib[idx][i] > 1.5]:
                    elem[i] += f"{int(lib[idx][i]):1d}"
        
        # Process Si
        hlpi = [lib[16][i] + lib[17][i] + lib[18][i] for i in range(length)]
        FIL = [i for i, x in enumerate(hlpi) if x > 0]
        if FIL:
            for i in FIL:
                Si[i] = 'Si'
            for i in [i for i in FIL if hlpi[i] > 1.5]:
                Si[i] += f"{int(hlpi[i]):1d}"

    # Lines 12256-12257: Process Hpos and eneg
    for i, x in enumerate(lib[9]):
        if x > 0:
            Hpos[i] = 'H+'
    for i, x in enumerate(lib[10]):
        if x > 0:
            eneg[i] = '+'

    # Line 12262: Combine all elements to form the formula
    form = [C13[i] + C12[i] + H1[i] + O18[i] + O17[i] + O16[i] + N15[i] + N14[i] + 
            Cl35[i] + Cl37[i] + F[i] + S32O4[i] + S34O4[i] + Si[i] + Hpos[i] + eneg[i] 
            for i in range(length)]

    return form

def get_eng_data(event: Any, file: str) -> Dict[str, Any]:
    # Line 12262-12263: Initialize info array
    info: List[str] = [''] * 9
    info[0] = 'H5 OPEN:'

    # Line 12265-12266: Construct file names
    file2 = file.replace('.h5', '_PTR.h5')
    
    # Line 12267-12269: Initialize error code and result structure
    err_code = 0
    s1 = {'duration': -9999}

    # Line 12273-12274: Open HDF5 file
    with h5py.File(file, 'r') as file_id:
        print(f'H5 FILE OPEN: {file}')

        # Line 12278-12297: Read datasets from file
        file_1_data = {'n_paths': get_par('#file1')}
        data_sets = [-1]
        for i in range(1, get_par('#file1') + 1):
            try:
                str_nam = f'H1_{i:02d}nam'
                str_dat = f'H1_{i:02d}dat'
                if str_nam in file_id:
                    file_1_data[str_nam] = file_id[str_nam][()]
                else:
                    file_1_data[str_nam] = '-9999'
                file_1_data[str_dat] = file_id[str_dat][()]
                data_sets.append(i)
            except:
                err_code = [err_code, i]

        # Line 12300-12304: Populate info array
        for i in range(1, len(data_sets)):
            str_dat = f'H1_{data_sets[i]:02d}dat'
            info[i] = str_dat

        # Line 12306-12307: Set instrument type
        instrument = 'TOF1000' if '/SPECdata/AverageSpec' in info else 'TOF8000'

        # Line 12311-12324: Read file attributes
        try:
            computer_id = 'nix'
            for attr_name, attr_value in file_id.attrs.items():
                if 'Computer ID' in attr_name:
                    computer_id = attr_value
        except:
            pass

        # Line 12325-12327: Check for VOCUS instrument
        if '2197a787-6c47-4c78-88d1-28fb6598bcc2' in computer_id or '3cef3dd0-f2f5-4216-95dd-9bdb829f125a' in computer_id:
            instrument = 'VOCUS'

        # Line 12333-12354: Read additional file if exists
        if os.path.exists(file2):
            with h5py.File(file2, 'r') as file_id2:
                file_2_data = {'n_paths': get_par('#file2')}
                data_sets2 = [-1]
                for i in range(1, get_par('#file2') + 1):
                    try:
                        str_nam = f'H2_{i:02d}nam'
                        str_dat = f'H2_{i:02d}dat'
                        if str_nam in file_id2:
                            file_2_data[str_nam] = file_id2[str_nam][()]
                        else:
                            file_2_data[str_nam] = '-9999'
                        file_2_data[str_dat] = file_id2[str_dat][()]
                        data_sets2.append(i)
                    except:
                        err_code = [err_code, i]

        # Line 12358-12445: Extract various parameters
        if max(data_sets) != -1:
            # Read file attributes
            creation_time = file_id.attrs.get('HDF5 File Creation Time', 'nix')
            
            if creation_time != 'nix':  # TOF1000
                cycle_time = float(file_id.attrs['Single Spec Duration (ms)']) / 1000  # in seconds
                sample_interval = max(float(file_id.attrs['Timebin width (ps)']) / 1e12)  # in seconds
                single_spectrum_time = float(file_id.attrs['Pulsing Period (ns)']) / 1000  # in microseconds
                pois_cor = 1
                pois_dead = 1
                extractions = max(cycle_time * 1e6 / single_spectrum_time)
                
                buf_times = file_1_data[f'H1_{data_sets[1]:02d}DAT'][3, :]
                buf_times = np.arange(len(buf_times)) * cycle_time + cycle_time
                sum_spec = file_1_data[f'H1_{data_sets[2]:02d}DAT'] * len(buf_times)
                computer_id = '-9999'
                single_ion_signal = 0.1
            else:  # TOF8000, QiTOF
                extractions = file_id['NbrWaveforms'][()]
                creation_time = file_id.attrs['HDF5_FILE_CREATION_TIME']
                sample_interval = file_id['fullspectra/sampleinterval'][()]
                pois_cor = file_id['RawData/HPTDC/PoissonCorrect'][()]
                computer_id = file_id.attrs['COMPUTER_ID']
                single_ion_signal = file_id['fullspectra/SINGLE_ION_SIGNAL'][()]
                pois_dead = file_id['RawData/HPTDC/PoissonTdcDeadTime'][()]
                
                buf_times = file_1_data[f'H1_{data_sets[1]:02d}DAT']
                sum_spec = file_1_data[f'H1_{data_sets[2]:02d}DAT']
                
                buf_times = buf_times.reshape(-1)
                time_step = np.median(np.diff(buf_times))
                buf_times = np.arange(len(buf_times)) * time_step

            # Ensure minimum signal threshold
            min_sig = 222 * max(buf_times)
            if not np.any(sum_spec > min_sig):
                raise ValueError("Signal too low")

            cycles = len(buf_times)
            buf_times = buf_times[:cycles]
            start_time = t09(creation_time)
            time_row = start_time + buf_times / 86400

        # Line 12449-12519: Process engineering data
        # Process engineering data
        filt_d = np.where((np.char.startswith(names_ptr, 'p_drift') | np.char.startswith(names_ptr, 'p-Drift')))[0]
        filt_u1 = np.where(np.char.startswith(names_ptr, 'Udrift'))[0]
        filt_u2 = np.where((np.char.startswith(names_ptr, 'Uq') | np.char.startswith(names_ptr, 'Udx')))[0]
        filt_t = np.where((np.char.startswith(names_ptr, 'Drift_Temperature') | np.char.startswith(names_ptr, 'T-Drift')))[0]

        cyc = np.zeros((1, cycles))
        pdrift = data[filt_d, :] if filt_d.size > 0 else cyc + get_par('p_drift_default')
        mist_ud = data[filt_u1, :][0] if filt_u1.size > 0 else cyc + get_par('u_drift_default')
        udrift = mist_ud + data[filt_u2, :] if filt_u2.size > 0 else mist_ud + get_par('udx_default')
        tdrift = data[filt_t, :] + 273.15 if filt_t.size > 0 else cyc + get_par('t_drift_default')

        pdrift = pdrift[0]
        udrift = udrift[0]
        tdrift = tdrift[0]

        if np.any(np.char.startswith(names_ptr, get_par('pdriftVOCUS'))):
            instrument = 'VOCUS'

        if instrument == 'VOCUS':
            pdrift = data[np.where(np.char.startswith(names_ptr, get_par('pdriftVOCUS')))[0][-1], :]
            udrift = data[np.where(np.char.startswith(names_ptr, get_par('udriftVOCUS')))[0][-1], :]
            tdrift = data[np.where(np.char.startswith(names_ptr, get_par('tdriftVOCUS')))[0][-1], :] + 273.15

        # Calculate parameters
        d = get_par('reactionlength')
        mu0 = get_par('reduced_mobility')
        mu = mu0 * (1013.25 / pdrift) * (tdrift / 273.15)
        trxn = d / ((mu * udrift) / d)
        n_molec = 24.63 * 298 * pdrift / (1013.25 * tdrift)
        e_over_n = 1e17 * (udrift / d) / (n_molec * 1e18)

        names = ['JulianDate', 'ReactionTime[s]', 'N [1E18 molec/cm3]', 'EoverN'] + names_ptr

        # Remove commas from names
        for i, name in enumerate(names):
            names[i] = name.replace(',', '_')

        # Adjust time_row and data lengths
        if time_row.shape[1] > trxn.shape[0]:
            time_row = time_row[:, :trxn.shape[0]]
        elif time_row.shape[1] < trxn.shape[0]:
            hlpp = np.zeros((1, trxn.shape[0]))
            hlpp[0, :time_row.shape[1]] = time_row
            time_row = hlpp

        data = np.vstack((time_row, trxn, n_molec, e_over_n, data))

        # Line 12523-12583: Create and return result structure
        s1 = {
            'names': names,
            'data': data.T,
            'FileName': file,
            'SumSpec': sum_spec,
            'cycles': cycles,
            'duration': max(buf_times),
            'pdrift': pdrift,
            'udrift': udrift,
            'tdrift': tdrift,
            'StartTime': start_time,
            'SampInt': sample_interval,
            'extractions': extractions,
            'PoisCor': pois_cor,
            'PoisDead': pois_dead,
            'instrument': instrument,
            'ComputerID': computer_id,
            'SingleIonSignal': single_ion_signal
        }

    return s1

def GetFileExport(event: Dict, destfolder: str, Name1: str) -> Dict[str, Union[List[float], List[List[float]], float, List[float]]]:
    # Lines 12589-12598: Get values from widgets
    massFL = event['Label_6']['uVALUE']
    massUL = event['Label_7']['uVALUE']
    cpsFL = event['base_1']['uVALUE']
    corrcpsFL = event['base_2']['uVALUE']
    cpsUL = event['base_3']['uVALUE']
    corrcpsUL = event['base_4']['uVALUE']
    ppbFL = event['base_5']['uVALUE']
    corrppbFL = event['base_6']['uVALUE']
    ppbUL = event['base_10']['uVALUE']
    corrppbUL = event['base_12']['uVALUE']

    # Line 12601: Get Filepar from widget
    Filepar = event['Label_5']['uVALUE']

    # Lines 12603-12625: Check if data needs to be loaded from files
    if max(max(cpsFL), max(corrcpsFL), max(cpsUL), max(corrcpsUL)) == -9999:
        # Load data from files
        massUL = read_csv(os.path.join(destfolder, f'Export/UL/IonList/MassIDs_{Name1}UL.csv'))
        massFL = read_csv(os.path.join(destfolder, f'Export/FL/IonList/MassIDs_{Name1}FL.csv'))
        corrcpsUL = read_float(os.path.join(destfolder, f'Export/UL/cps/ocorr/corrcps{Name1}UL.fdt'))
        cpsUL = read_float(os.path.join(destfolder, f'Export/UL/cps/cps{Name1}UL.fdt'))
        corrcpsFL = read_float(os.path.join(destfolder, f'Export/FL/cps/ocorr/corrcps{Name1}FL.fdt'))
        cpsFL = read_float(os.path.join(destfolder, f'Export/FL/cps/cps{Name1}FL.fdt'))
        corrppbUL = read_float(os.path.join(destfolder, f'Export/UL/ppb/ocorr/corrppb{Name1}UL.fdt'))
        ppbUL = read_float(os.path.join(destfolder, f'Export/UL/ppb/ppb{Name1}UL.fdt'))
        corrppbFL = read_float(os.path.join(destfolder, f'Export/FL/ppb/ocorr/corrppb{Name1}FL.fdt'))
        ppbFL = read_float(os.path.join(destfolder, f'Export/FL/ppb/ppb{Name1}FL.fdt'))

        # Update widget values
        event['Label_6']['uVALUE'] = massFL
        event['Label_7']['uVALUE'] = massUL
        event['base_1']['uVALUE'] = cpsFL
        event['base_2']['uVALUE'] = corrcpsFL
        event['base_3']['uVALUE'] = cpsUL
        event['base_4']['uVALUE'] = corrcpsUL
        event['base_5']['uVALUE'] = ppbFL
        event['base_6']['uVALUE'] = corrppbFL
        event['base_10']['uVALUE'] = ppbUL
        event['base_12']['uVALUE'] = corrppbUL

    # Lines 12627-12628: Return a dictionary with all the data
    return {
        'massUL': massUL,
        'massFL': massFL,
        'corrcpsUL': corrcpsUL,
        'cpsUL': cpsUL,
        'corrcpsFL': corrcpsFL,
        'cpsFL': cpsFL,
        'corrppbUL': corrppbUL,
        'ppbUL': ppbUL,
        'corrppbFL': corrppbFL,
        'ppbFL': ppbFL,
        'resolution': Filepar['resolution'],
        'peakshape': Filepar['peakshape']
    }

def getPrimIons(masses: np.ndarray, cps: np.ndarray, mode: int) -> Dict[str, np.ndarray]:
    # Line 12632: Function definition
    
    # Lines 12633-12650: H3O+ mode
    if mode == 0:
        # H3O+   m21factor = 1/0.206% = 487 m38factor = 99.926/0.074 = 1350 , m39factor = 99.588/0.412 = 242
        Prim1 = cps[np.abs(masses - 21.0221).argmin(), :]
        if np.median(Prim1) > 100:
            Prim1 *= 487
        else:
            Prim1 = cps[np.abs(masses - 19.018).argmin(), :]
            print('m19 used (because m21 below 100 cps)')
        
        if getpar('m38') == 1:
            Prim2 = cps[np.abs(masses - 38.0326).argmin(), :] * 645
        else:
            Prim2 = cps[np.abs(masses - 39.0327).argmin(), :] * 242
        
        if np.min(np.abs(masses - 38.0326)) > 0.005:
            Prim2 = np.zeros_like(Prim2)
            Prim2 = cps[np.abs(masses - 39.0327).argmin(), :] * 242
            print('Warning: m38.033 not detected')
        
        Prim1 = Prim1[0, :]
        Prim2 = Prim2[0, :]
    
    # Lines 12651-12656: NO+ mode
    elif mode == 5:
        Prim1 = cps[np.abs(masses - 30.994).argmin(), :] * 269
        Prim2 = cps[np.abs(masses - 47.9966).argmin(), :] * 242
        Prim1 = Prim1[0, :]
        Prim2 = Prim2[0, :]
    
    # Lines 12657-12662: O2+ mode
    elif mode == 9:
        Prim1 = cps[np.abs(masses - 33.9935).argmin(), :] * 242
        Prim2 = np.zeros_like(Prim1)
        Prim1 = Prim1[0, :]
        Prim2 = Prim2[0, :]
    
    # Lines 12664-12665: Check average primary ions threshold
    if np.mean(Prim1) < getpar('AveragePrimIonsThresh'):
        Prim1 = np.ones_like(Prim1) * np.mean(Prim1)
    
    # Lines 12667-12668: Return result as a dictionary
    return {'A': Prim1, 'B': Prim2}

def getpar(par_name: str) -> Union[float, str]:
    # Line 12673: Construct the parameter name pattern
    par_name2 = f'_{par_name}='
    
    # Lines 12674-12677: Read the parfile path
    with open('C:/PTRwid/parfile.txt', 'r') as lun:
        path = lun.readline().strip()
    
    # Line 12678: Read the parameter file content
    with open(path, 'r') as f:
        ptrwidpar = f.readlines()
    
    # Line 12679: Find the row index of the parameter
    row_index = next((i for i, line in enumerate(ptrwidpar) if par_name2 in line), -1)
    
    if row_index == -1:
        return None  # Parameter not found
    
    # Lines 12680-12681: Extract the value part of the parameter
    line = ptrwidpar[row_index]
    value_start = line.index('=') + 1
    value_str = line[value_start:].strip()
    
    # Line 12682: Convert to float if possible, otherwise return as string
    try:
        return float(value_str)
    except ValueError:
        return value_str

def getstr(par_name: str) -> str:
    # Line 12687: Construct parameter name with underscore prefix and equals sign
    par_name2 = f'_{par_name}='

    # Lines 12688-12691: Read the parfile path from C:/PTRwid/parfile.txt
    with open('C:/PTRwid/parfile.txt', 'r') as lun:
        path = lun.readline().strip()

    # Line 12692: Read the parameter file content
    with open(path, 'r') as f:
        ptrwidpar = f.readlines()

    # Line 12693: Find the row index where the parameter name matches
    row_index = next(i for i, line in enumerate(ptrwidpar) if par_name2 in line)

    # Lines 12694-12695: Find the position of the equals sign in the matching line
    line_index = ptrwidpar[row_index].index('=') + 1

    # Line 12696: Extract the parameter value
    str_value = ptrwidpar[row_index][line_index:].strip()

    # Line 12697: Return the trimmed string value
    return str_value

def get_inst_id() -> str:
    """
    Read and return the instrument ID from a file.
    Corresponds to lines 12700-12705 in the original IDL code.
    """
    inst_id = ''
    with open('C:/PTRwid/InstrumentID.txt', 'r') as file:
        inst_id = file.read().strip()
    return inst_id

def get_tr(unique_campaign_name: str, zeit: float = 0) -> List[float]:
    """
    Retrieve transmission parameters for a given campaign and time.
    Corresponds to lines 12707-12733 in the original IDL code.
    """
    with open('C:/PTRwid/TransmissionArchive.txt', 'r') as file:
        tr_arch = file.readlines()
    
    par_name = f'_{unique_campaign_name}='
    row_index = next((i for i, line in enumerate(tr_arch) if par_name in line), -1)
    
    if row_index == -1:
        row_index = next(i for i, line in enumerate(tr_arch) if '_default=' in line)
        print('WARNING: UniqueCampaignName not found. Using default parameter')
    
    line = tr_arch[row_index]
    str_values = line.split('=')[1].strip().split(',')
    sets = [list(map(float, values.split(','))) for values in str_values]
    
    filt = [i for i, set_values in enumerate(sets) if set_values[0] < zeit]
    if filt:
        set_values = sets[max(filt)]
    else:
        set_values = sets[0]
    
    if abs(zeit - set_values[0]) > 30:
        print('WARNING: nearest transmission more than 30 days away')
    
    return [value for value in set_values[1:] if value != -1]

def add_tr(unique_campaign_name: str, set_values: List[float]) -> None:
    """
    Add new transmission parameters to the archive.
    Corresponds to lines 12736-12749 in the original IDL code.
    """
    with open('C:/PTRwid/TransmissionArchive.txt', 'r') as file:
        tr_arch = file.readlines()
    
    new_line = f'_{unique_campaign_name}=' + ','.join(map(str, set_values)) + '\n'
    tr_arch.append(new_line)
    
    with open('C:/PTRwid/TransmissionArchive.txt', 'w') as file:
        file.writelines(tr_arch)

def isodist(lib: List[float]) -> Dict[str, List[float]]:
    # Lines 12752-12768: Define isotope probabilities and masses
    C12p, C13p = 0.989, 0.011
    H1p = 1
    O16p, O17p, O18p = 0.9976, 0.0004, 0.0020
    N14p, N15p = 0.9963, 0.0037
    H_posp, e_negp = 1, 1
    Cl35p, Cl37p = 0.7577, 0.2423
    Fp = 1
    S32p, S33p, S34p = 0.9504, 0.0075, 0.0421

    C12, C13 = 12, 13.003355
    H1 = 1.007825
    O16, O17, O18 = 15.994915, 16.999131, 17.99916
    N14, N15 = 14.003074, 15.000108
    H_pos, e_neg = 1.007276467, -0.000548533
    Cl35, Cl37 = 34.968852, 36.965903
    F = 18.998403
    S32O4 = 31.97207 + 4 * O16
    S33O4 = 32.971456 + 4 * O16
    S34O4 = 33.967866 + 4 * O16

    # Lines 12771-12773: Calculate number of atoms
    nC = lib[0] + lib[1]
    nN = lib[6] + lib[7]
    nCl = lib[10] + lib[11]

    # Lines 12775-12777: Calculate probabilities
    P13C = nC * C13p * C12p**(nC-1)
    P15N = nN * N15p * N14p**(nN-1)
    P37Cl = nCl * Cl37p * Cl35p**(nCl-1)

    # Lines 12780-12781: Initialize result arrays
    masses: List[float] = [0]
    disp: List[float] = [0]

    # Lines 12782-12799: Main calculation loop
    for c in range(2):
        for nnn in range(2):
            for cl in range(2):
                if c <= nC and nnn <= nN and cl <= nCl:
                    n = [nC-c, c] + lib[2:6] + [nN-nnn, nnn] + lib[8:10] + [nCl-cl, cl] + lib[12:15]
                    entry = (C12*n[0] + C13*n[1] + H1*n[2] + O16*n[3] + O17*n[4] + O18*n[5] +
                             N14*n[6] + N15*n[7] + Cl35*n[10] + Cl37*n[11] + F*n[12] +
                             S32O4*n[13] + S34O4*n[14] + H_pos*n[8] + e_neg*n[9])
                    PC = P13C if c == 1 else (1 - P13C) if nC >= 1 else 1
                    PN = P15N if nnn == 1 else (1 - P15N) if nN >= 1 else 1
                    PCl = P37Cl if cl == 1 else (1 - P37Cl) if nCl >= 1 else 1
                    Ptot = PC * PN * PCl
                    disp.append(Ptot)
                    masses.append(entry)

    # Lines 12800-12803: Filter and sort results
    filtered = [(m, d) for m, d in zip(masses, disp) if d >= 0.002]
    filtered.sort(key=lambda x: x[0])
    masses, disp = zip(*filtered)

    # Line 12808: Create and return result structure
    return {'masses': list(masses), 'dist': list(disp)}

def LoadMassRange(file: str, MassLow: float, MassHigh: float, a: float, t0: float, ex: float, SampInt: float) -> Union[List[float], int]:
    # Line 12825-12827: Initialize variables
    ErrCode: int = 0
    ptrData: Union[List[float], int] = -9999
    locatie: str = '/FullSpectra/TofData'

    # Lines 12828-12837: Error handling
    try:
        # Line 12839-12844: Open HDF5 file and dataset
        with h5py.File(file, 'r') as fileid:
            dataset_id1 = fileid[locatie]
            dataspace_id1 = dataset_id1.id.get_space()
            Dims: Tuple[int, ...] = dataspace_id1.get_simple_extent_dims()

            # Line 12846-12848: Calculate start and width
            Start: int = max(1, int(m2t(MassLow, a, t0, ex, SampInt, mass=True)))
            Width: int = int(m2t(MassHigh, a, t0, ex, SampInt, mass=True) - Start)

            # Lines 12850-12859: Set up hyperslab selection
            if len(Dims) >= 2:
                ddd: List[int] = [0] * (len(Dims) - 1)
                eee: List[int] = list(Dims[1:])
                start2: List[int] = [Start] + ddd
                count: List[int] = [Width] + eee
            else:
                start2: List[int] = [Start]
                count: List[int] = [Width]

            # Lines 12860-12862: Select hyperslab and read data
            dataset_id1.select(start2, count)

            ptrData = dataset_id1[tuple(slice(start2[i], start2[i]+count[i]) for i in range(len(start2)))]

    except Exception as e:
        # Lines 12830-12836: Handle errors
        if locatie == '/FullSpectra/TofData':
            ErrCode = 1
            locatie = '/SPECdata/Intensities'
            if ErrCode == 1:
                print(f'Error: {str(e)}')
        else:
            ptrData = -9999
            print('HDF5 FILE CORRUPTED 3')

    # Line 12871: Return the data
    return ptrData

def getSplit(file: str) -> dict:
    # Line 12875-12877: Initialize variables
    ErrCode: int = 0
    ptrData: float = -9999
    locatie: str = '/FullSpectra/TofData'

    # Lines 12878-12888: Error handling (simplified for Python)
    try:
        with h5py.File(file, 'r') as fileid:
            dataset = fileid[locatie]
    except KeyError:
        locatie = '/SPECdata/Intensities'
        try:
            with h5py.File(file, 'r') as fileid:
                dataset = fileid[locatie]
        except KeyError:
            print(f"Error: Unable to open dataset at {locatie}")
            return {'split': 1, 'start2': [], 'count': [], 'locatie': locatie}

    # Lines 12891-12896: Get dataset dimensions
    Dims: Tuple[int, ...] = dataset.shape
    n_dims: int = len(Dims)
    JunkSize: float = getpar('JunkSize')  # Assuming getpar function is defined elsewhere

    # Lines 12898-12900: Calculate total points
    timebins: int = Dims[0]
    pnts: int = timebins
    for i in range(1, n_dims):
        pnts *= Dims[i]

    # Lines 12902-12938: Calculate split if necessary
    if pnts > JunkSize:
        split: int = 0
        counter: int = 0
        while split == 0:
            pntsSUB: int = timebins
            if n_dims > 2.5:
                for i in range(2, n_dims - 1 - counter):
                    pntsSUB *= Dims[i]
            if pntsSUB < JunkSize:
                max_split: int = Dims[n_dims - 1 - counter]
                spliti: int = int(np.ceil(pntsSUB * max_split / JunkSize))
                counti: int = int(round(max_split / spliti))

                start2: List[List[int]] = [list(Dims) for _ in range(spliti)]
                count: List[List[int]] = [list(Dims) for _ in range(spliti)]
                for i in range(spliti):
                    start2[i][n_dims - 1 - counter] = i * counti
                    if i == spliti - 1:
                        count[i][n_dims - 1 - counter] = Dims[n_dims - 1 - counter] - (spliti - 1) * counti
                    else:
                        count[i][n_dims - 1 - counter] = counti

                split = spliti
            else:
                counter += 1
    else:
        split: int = 1
        start2: List[List[int]] = [list(Dims)]
        count: List[List[int]] = [list(Dims)]

    # Lines 12946-12950: Return results
    return {
        'split': split,
        'start2': start2,
        'count': count,
        'locatie': locatie
    }

def LoadCycles(file: str, locatie: str, strt: List[int], cnt: List[int]) -> Tuple[int, List[float]]:
    """
    Load cycles from an HDF5 file.
    
    Args:
    file (str): Path to the HDF5 file
    locatie (str): Dataset location in the HDF5 file
    strt (List[int]): Start indices for hyperslab selection
    cnt (List[int]): Count values for hyperslab selection
    
    Returns:
    Tuple[int, List[float]]: Error code and loaded data
    """
    # Line 12957-12958: Initialize error code and data
    ErrCode = 0
    ptrData = -9999

    # Line 12959-12960: Check for errors
    if ErrCode == 0:
        try:
            # Lines 12962-12965: Open HDF5 file and dataset
            with h5py.File(file, 'r') as fileid:
                dataset = fileid[locatie]
                
                # Lines 12966-12967: Get dataset dimensions
                Dims = dataset.shape
                
                # Lines 12974-12976: Select hyperslab and read data
                data = dataset[tuple(slice(s, s+c) for s, c in zip(strt, cnt))]
                
                # Line 12977: Rearrange data (implement rearrange function as needed)
                ptrData = rearrange(data)

        except Exception as e:
            print(f"Error reading HDF5 file: {e}")
            ErrCode = 1
            ptrData = [-9999]

    else:
        # Lines 12991-12993: Handle error case
        ptrData = [-9999]
        print('HDF5 FILE CORRUPTED 3')

    # Line 12995: Return results
    return ErrCode, ptrData

def m2t(value: Union[float, np.ndarray], a: float, t0: float, ex: float, SampInt: float, 
        time: Optional[bool] = False, mass: Optional[bool] = False) -> Union[float, np.ndarray]:
    """
    Convert TOF time scale to TOF mass scale and vice versa.
    
    # Lines 13001-13002: Function definition and docstring
    """
    
    # Lines 13004-13005: Check if SampInt is defined
    if 'SampInt' not in locals():
        SampInt = 1e-10

    # Lines 13010-13016: Determine conversion direction based on keywords
    if time:
        # Force time to mass conversion
        value2 = np.exp(np.log((value * SampInt / 1e-10 - t0) / a) / ex)
    elif mass:
        # Force mass to time conversion
        value2 = (a * value**ex + t0) * 1e-10 / SampInt
    else:
        # Automatic determination based on value range
        if np.min(value) > 3500:
            # Time to mass conversion
            value2 = np.exp(np.log((value * SampInt / 1e-10 - t0) / a) / ex)
        else:
            # Mass to time conversion
            value2 = (a * value**ex + t0) * 1e-10 / SampInt

    # Line 13023: Return the converted value
    return value2

def MakeCsv(name: str, data: Union[List[List[str]], List[str]]) -> None:
    """
    Write data to a CSV file.
    
    Args:
    name (str): The filename to write to.
    data (Union[List[List[str]], List[str]]): The data to write to the file.
    """
    # Lines 13028-13036: Error handling
    try:
        # Lines 13038-13052: Main CSV writing logic
        with open(name, 'w', newline='') as f:
            
            
            # Determine if data is 2D or 1D
            if isinstance(data[0], list):
                # 2D data
                for row in data:
                    # Convert all elements to strings and remove leading/trailing whitespace
                    row = [str(item).strip() for item in row]
                    # Add commas between columns, except for the last column
                    for i in range(len(row) - 1):
                        row[i] += ','
                    f.write(f"{row}\n")
            else:
                # 1D data
                # Convert all elements to strings and remove leading/trailing whitespace
                data = [str(item).strip() for item in data]
                # Add commas between items, except for the last item
                for i in range(len(data) - 1):
                    data[i] += ','
                f.write(f"{row}\n")
        
        print(f"CSV file '{name}' has been created successfully.")
    
    except Exception as e:
        # Lines 13030-13035: Error reporting
        print(f"Error in MakeCsv: {str(e)}")
        print("MAKECSV not possible (disabled) :-((")

def read_csv(file: str, semicolon: bool = False, decimal_comma: bool = False, return_string: bool = False) -> Union[List[List[float]], List[List[str]]]:
    """
    Read a CSV file and return its contents as a list of lists.

    Args:
    file (str): Path to the CSV file
    semicolon (bool): If True, use semicolon as separator instead of comma
    decimal_comma (bool): If True, replace commas with dots in numeric values
    return_string (bool): If True, return data as strings instead of floats

    Returns:
    List[List[Union[float, str]]]: The contents of the CSV file
    """
    # Line 25468: Set separator based on semicolon flag
    separator = ';' if semicolon else ','

    # Lines 25470-25471: Check if file exists
    if not os.path.exists(file):
        return [[-9999]]  # Return default value if file doesn't exist

    # Lines 25472-25477: Read file contents
    with open(file, 'r') as f:
        data = f.readlines()

    # Lines 25479-25480: Determine number of columns
    cols = len(data[-1].split(separator))

    # Lines 25481-25482: Initialize result list
    result = []

    # Lines 25483-25495: Process each line
    for line in data:
        # Replace decimal comma if needed
        if decimal_comma:
            line = line.replace(',', '.')
        
        # Split line and process values
        values = line.strip().split(separator)
        if len(values) == cols:
            if return_string:
                result.append(values)
            else:
                try:
                    result.append([float(v) for v in values])
                except ValueError:
                    # Skip lines that can't be converted to float
                    continue

    return result

def read_csv_str(file: str, long: bool = False, semicolon: bool = False, variable: bool = False) -> List[Union[str, List[str]]]:
    # Line 25498-25499: Set separator based on semicolon flag
    separator = ';' if semicolon else ','

    # Line 25501-25502: Check if file exists
    if os.path.exists(file):
        # Line 25503-25506: Open file and read lines
        with open(file, 'r') as f:
            data = f.readlines()
        
        # Line 25507-25508: Strip newlines and use appropriate format
        format_str = '(A200)' if long else '(A30)'
        data = [line.strip() for line in data]

        # Line 25510-25511: Determine number of columns
        cols = max(len(line.split(separator)) for line in data)

        if variable:
            # Line 25514-25519: Handle variable-length lines
            data2 = [[''] * cols for _ in range(len(data))]
            for i, line in enumerate(data):
                parts = line.split(separator)
                data2[i][:len(parts)] = parts
        else:
            # Line 25521-25524: Handle fixed-length lines
            data2 = [line.split(separator) for line in data]
            if len(data) > 1:
                data2 = list(map(list, zip(*data2)))  # Transpose data

    else:
        # Line 25529-25531: Return error value if file doesn't exist
        data2 = ['-9999']

    return data2

def read_index_file(index_file: str) -> Dict[str, Union[List[str], str]]:
    # Line 25540-25541: Check if file exists and has correct extension
    if os.path.exists(index_file) and index_file.endswith('.ind'):
        # Line 25542-25546: Read the entire file
        with open(index_file, 'r') as file:
            index_par = file.readlines()
        
        # Line 25548-25550: Extract interval
        interval = next((line.split('=')[1].strip() for line in index_par if '_Textt0=' in line), '')
        
        # Line 25553-25557: Initialize arrays
        vector_ind = [''] * 16
        cond_ind = [''] * 16
        value_ind = [''] * 16
        if_ind = [''] * 16
        else_ind = [''] * 16
        
        # Line 25559-25575: Extract values for each index
        for i in range(1, 17):
            ii = str(i) if i < 10 else f"{i}"
            vector_ind[i-1] = next((line.split('=')[1].strip() for line in index_par if f'_Tex{ii}=' in line), '')
            cond_ind[i-1] = next((line.split('=')[1].strip() for line in index_par if f'_Dropp{ii}=' in line), '')
            value_ind[i-1] = next((line.split('=')[1].strip() for line in index_par if f'_Text{ii}=' in line), '')
            if_ind[i-1] = next((line.split('=')[1].strip() for line in index_par if f'_Textt{ii}=' in line), '')
            else_ind[i-1] = next((line.split('=')[1].strip() for line in index_par if f'_Texttt{ii}=' in line), '')
        
        # Line 25577-25586: Extract additional values
        ind_value = [''] * 4
        values_2_add = [''] * 4
        for i in range(15, 19):
            ii = f"{i}"
            ind_value[i-15] = next((line.split('=')[1].strip() for line in index_par if f'_Dext{ii}=' in line), '')
            values_2_add[i-15] = next((line.split('=')[1].strip() for line in index_par if f'_Dextt{ii}=' in line), '')
        
        # Line 25589-25590: Create and return result structure
        return {
            'VectorInd': vector_ind,
            'CondInd': cond_ind,
            'ValueInd': value_ind,
            'ifInd': if_ind,
            'elseInd': else_ind,
            'IndValue': ind_value,
            'Values2add': values_2_add,
            'interval': interval
        }
    else:
        # Line 25594: Return default structure if file doesn't exist or has wrong extension
        return {'VectorInd': -9999}

def ReadFloat(file: str) -> Union[List[float], float]:
    """
    Read float data from a file.
    
    Lines 25603-25625: Main function logic
    """
    if os.path.exists(file):
        with open(file, 'r') as f:
            data = f.readlines()
        
        # Lines 25612-25613: Replace -NaN values
        data = ['-9.999000E+003' if line.strip() == '-NaN' else line for line in data]
        data = ['-9.999000000E+003' if line.strip() == '-NaN' else line for line in data]
        
        # Lines 25614-25617: Convert to float and reshape
        data2 = [float(line) for line in data]
        n_dim = int(data2[0])
        if n_dim > 0.5:
            data2 = [data2[n_dim+1:]]
        else:
            data2 = data2[2]
    else:
        data2 = -9999
    
    return data2

def rearrange(Data: List[float]) -> List[List[float]]:
    """
    Rearrange data into a 2D list.
    
    Lines 25629-25634: Function logic
    """
    mist = [len(Data)]  # Assuming Data is a 1D list
    n_dims = 1  # Assuming 1D input
    mist3 = 1
    for i in range(1, n_dims):
        mist3 *= mist[i]
    Data = [Data[:mist[0]]] * mist3
    return Data

def starttimes(files: List[str], ID: str) -> Union[List[float], float]:
    """
    Calculate start times for given files based on their filenames.
    
    Args:
    files (List[str]): List of file names
    ID (str): Identifier string in the file names
    
    Returns:
    Union[List[float], float]: List of start times or 0 if no valid times found
    """
    # Line 25659: Get dimensions of files
    mist = len(files)
    
    # Line 25660: Get length of ID
    n = len(ID)
    
    # Lines 25661-25662: Transpose files if necessary
    if mist > 2:
        files = [file for file in zip(*files)]
    
    # Lines 25663-25668: Extract date and time components from filenames
    Year = [file[file.index(ID)+n-2:file.index(ID)+n+2] for file in files]
    Month = [file[file.index(ID)+n+3:file.index(ID)+n+5] for file in files]
    Day = [file[file.index(ID)+n+6:file.index(ID)+n+8] for file in files]
    Hour = [file[file.index(ID)+n+9:file.index(ID)+n+11] for file in files]
    Minute = [file[file.index(ID)+n+12:file.index(ID)+n+14] for file in files]
    Second = [file[file.index(ID)+n+15:file.index(ID)+n+17] for file in files]
    
    # Lines 25669-25670: Calculate t09 (days since Jan 1, 2009)
    if max(Year) < 0.1:
        t09 = 0
    else:
        t09 = [(datetime(int(y), int(m), int(d), int(h), int(mi), int(s)) - 
                datetime(2009, 1, 1)).total_seconds() / (24 * 3600) 
               for y, m, d, h, mi, s in zip(Year, Month, Day, Hour, Minute, Second)]
    
    return t09

def str2vec(str_input: str) -> Union[List[float], float]:
    """
    Convert a string of comma-separated values to a list of floats.
    
    Args:
    str_input (str): Input string
    
    Returns:
    Union[List[float], float]: List of float values or -9999 if invalid input
    """
    # Lines 25675-25676: Check if input is a valid number
    NoNumber = any(char.isalpha() for char in str_input.lower())
    
    # Lines 25679-25681: Convert string to vector of floats
    if str_input == '':
        vec = -9999
    elif NoNumber:
        vec = -9999
    else:
        vec = [float(x) for x in str_input.split(',')]
    
    return vec

def t09(CreationTime: str, format6: bool = False) -> float:
    """
    Convert various date-time formats to days since January 1, 2009.
    
    Args:
    CreationTime (str): The input date-time string.
    format6 (bool): Flag to use format6 parsing.
    
    Returns:
    float: Days since January 1, 2009.
    """
    # Reference date
    ref_date = datetime(2009, 1, 1)
    
    # Parse different formats
    if '/' in CreationTime and 'h' not in CreationTime:
        # Format 1 and 1b
        if 'PM' in CreationTime or 'AM' in CreationTime:
            dt = datetime.strptime(CreationTime, "%m/%d/%Y, %I:%M:%S %p")
        else:
            try:
                dt = datetime.strptime(CreationTime, "%d/%m/%Y %H:%M:%S")
            except ValueError:
                dt = datetime.strptime(CreationTime, "%m/%d/%Y %H:%M:%S")
    elif '-' in CreationTime and not format6:
        # Format 2
        dt = datetime.strptime(CreationTime, "%d-%b-%y %I:%M:%S %p")
    elif '/' in CreationTime and 'h' in CreationTime:
        # Format 3
        dt = datetime.strptime(CreationTime, "%d/%m/%Y %Hh %Mm %Ss")
    elif '.' in CreationTime and len(CreationTime.split('.')[0]) == 4:
        # Format 5
        dt = datetime.strptime(CreationTime, "%Y.%m.%d-%Hh%Mm%Ss")
    elif format6:
        # Format 6
        dt = datetime.strptime(CreationTime, "%Y-%m-%d %H:%M:%S")
    else:
        # Format 4
        dt = datetime.strptime(CreationTime, "%d.%m.%Y %H:%M:%S")
    
    # Calculate days since reference date
    days_since_ref = (dt - ref_date).total_seconds() / (24 * 3600)
    return days_since_ref

def t09str(date: Union[float, List[int]]) -> str:
    """
    Convert days since January 1, 2009 or a date list to a formatted string.
    
    Args:
    date (Union[float, List[int]]): Days since Jan 1, 2009 or [year, month, day, hour, minute, second].
    
    Returns:
    str: Formatted date string.
    """
    ref_date = datetime(2009, 1, 1)
    
    if isinstance(date, (int, float)):
        dt = ref_date + timedelta(days=date)
    else:
        dt = datetime(*date)
    
    return dt.strftime("%Y.%m.%d-%Hh%Mm%Ss")

def kml(compound_name: str, startTIME: str, endTIME: str) -> None:
    # Read CSV data
    data = read_csv('D:/Pamarcmip/polar 5/test13.csv', return_string=True)

    # Filter data based on time range
    filt = [i for i, row in enumerate(data) if float(startTIME) < float(row[0]) < float(endTIME)]
    data2 = [[float(val) for val in data[i]] for i in filt]

    # Extract longitude, latitude, altitude, and VMR data
    lon: List[float] = [row[1] for row in data2]
    lat: List[float] = [row[2] for row in data2]
    alt: List[float] = [row[3] for row in data2]
    VMR: List[float] = [row[data[0].index(compound_name)] for row in data2]

    # Define color schemes
    cols: List[str] = ['64000000', '64060055', '640D00AA', '641300F5', '641461FF', '641485FF', '6414BAFF', '6414E7FF', '6414F0B3', '6414F009', '64A7F9A0']
    cols = ['64000000', '64580007', '6485000B', '64BA000F', '64E70013', '649E2C6B', '645D50B0', '643466DC', '641C73F6', '641435FF', '641404FF']

    # Calculate quantiles for color mapping
    sep: List[float] = [quantile(VMR, 0.05 + i/11.0) for i in range(11)]

    # Initialize KML file structure
    dims: Tuple[int, int] = (len(data2), len(data2[0]))
    klmfile: List[str] = [''] * (8 + 14 * dims[0])
    klmfile[0] = '<?xml version="1.0" encoding="UTF-8"?>'
    klmfile[1] = '<kml xmlns="http://www.opengis.net/kml/2.2">'
    klmfile[2] = '  <Document>'
    klmfile[3] = ' <Folder>'
    klmfile[4] = '      <name>Trail</name>'
    klmfile[5 + 14 * dims[0]] = '    </Folder>'
    klmfile[6 + 14 * dims[0]] = '  </Document>'
    klmfile[7 + 14 * dims[0]] = '</kml>'

    # Generate KML content for each data point
    for i in range(dims[0] - 1):
        colind = next((j for j, s in enumerate(sep) if VMR[i] <= s), 10)
        colstr = cols[colind]
        coostr = f"{lon[i]:.6f},{lat[i]:.6f},{10.0*alt[i]:.3f} {lon[i+1]:.6f},{lat[i+1]:.6f},{10.0*alt[i+1]:.3f}"
        
        klmfile[5 + i*14:19 + i*14] = [
            '  <Placemark>',
            '        <Style>',
            '          <LineStyle>',
            f'            <color>{colstr}</color>',
            '            <width>3</width>',
            '          </LineStyle>',
            '       </Style>',
            '        <MultiGeometry>',
            '          <LineString>',
            '            <altitudeMode>absolute</altitudeMode>',
            f'            <coordinates>{coostr}</coordinates>',
            '          </LineString>',
            '        </MultiGeometry>',
            '      </Placemark>'
        ]

    # Print KML content
    for line in klmfile:
        print(line)
    print('haha')



