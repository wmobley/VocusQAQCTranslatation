import vocustools.auxilary_routines as aux
import os
import vocustools.unified_mass_list as uml
aux.setup_ptrwid_files('v003', date='jul_01_2021')

aux.create_ptrwid_file_structure('idl_example/Test_Files')




h5_files = []
for root, dirs, files in os.walk('idl_example/Test'):
    for file in files:
        if file.endswith('.h5'):
            relative_path = os.path.relpath(os.path.join(root, file), './')
            h5_files.append(relative_path)

print(h5_files)
names=['sig19','sig37','sig19TR','sig37TR','R37/17TR',    'R32/19TR','R78/79TR','R120/121TR','R106/107TR','Tdrift',    'pdrift','Udrift','Udx','Efield','treact',    'EoverN','','','','', 
        'nS_33','nS_42','nS_45','nS_59','nS_87',    'nS_69','nS_41','nS_71','nS_73','nS_79',    'nS_107','nS_121','nS_137','nS_81','nS_133',    'nS_181','nS_223','nS_207','nS_297','nS_281',    'nS_371','nS_355',
        'n59S_33','n59S_42','n59S_45','n59S_59','n59S_87',    'n59S_69','n59S_41','n59S_71','n59S_73','n59S_79',    'n59S_107','n59S_121','n59S_137','n59S_81','n59S_133',    'n59S_181','n59S_223','n59S_207','n59S_297','n59S_281',    'n59S_371','n59S_355',
        'rawS_33','rawS_42','rawS_45','rawS_59','rawS_87',    'rawS_69','rawS_41','rawS_71','rawS_73','rawS_79',    'rawS_107','rawS_121','rawS_137','rawS_81','rawS_133',    'rawS_181','rawS_223','rawS_207','rawS_297','rawS_281',    'rawS_371','rawS_355']

uml.uni_mass_list(h5_files, names, 'idl_example/Test_Files/MASSLIST')