''' 
Valentin Gonay
19/01/24
'''



import os
import pkg_resources
import importlib.util

GREEN = '\033[92m'
ORANGE = '\033[93m'
RED = '\033[91m'
RESET = '\033[0m'


def check_library(library: dict):
    '''Check if a library can be imported and optionally if it matches a specific version
    
    :param library: The library that will be used for the import. Need to have 'dep_name' (the \
        name of the installed library), 'import_name' (the name use to import the library in a \
        python script), and can have a 'version' (library version)
    :type library: dict 
    '''

    # Check if the module can be imported (standard library or installed package)
    dependency_name = library['dep_name']
    import_name = library['import_name']
    version = library['version']
    try:
        spec = importlib.util.find_spec(import_name)
        if spec is None:
            raise ImportError
        module_found = True
    except ImportError:
        module_found = False

    if module_found:
        # If version check is not required, print success message
        if not version:
            print(f"{dependency_name}"+GREEN+" is installed."+RESET)
            return
    else:
        print(f"{dependency_name}"+RED+" is not installed."+RESET)
        return

    # Check for the version using pkg_resources for installed packages
    try:
        dist = pkg_resources.get_distribution(dependency_name)
        if dist.version == version:
            print(f"{dependency_name} {version}"+GREEN+" is installed."+RESET)
        else:
            print(f"{dependency_name} {version}"+ORANGE+" is not in a tested version."+RESET)
            print(f"Installed version is: {dist.version}")
    except pkg_resources.DistributionNotFound:
        print(f"{dependency_name}"+RED+" is not installed."+RESET)
    except pkg_resources.VersionConflict:
        print(f"{dependency_name} has a version conflict."+RED+" Please check your",
              "installations."+RESET)
    except Exception as e:
        print(f"An error occurred: {e}")



def check_install():
    '''Check installation of all folders, files and library needed to run a prediction
    
    :return: True if there is a problem in the installation, else False
    :rtype: bool
    '''

    is_pb = False
    current_dir = os.path.dirname(__file__).replace('/utils','')

    # Check files and folders integrity
    folder_list = [
        'data', 
        'results', 
        'utils', 
        'utils/iupred3', 
        'utils/iupred3/.idea', 
        'utils/iupred3/data'
        ]
    for folder_name in folder_list:
        folder_path = os.path.join(current_dir,folder_name)
        if not os.path.exists(folder_path):
            print(f"Folder {folder_path}" +RED+ "\tis missing..."+RESET)
            is_pb = True
        else:
            print(f"folder {folder_path}"+GREEN+ "\tOK"+RESET)

    files_list = [
        'CB_RF_pred.py',
        'utils/progress_bar.py', 'utils/fold_pred.py', 
        'utils/get_features.py', 'utils/progress_bar.py',
        'utils/iupred3/iupred3_lib.py', 'utils/iupred3/iupred3.py',
        'data/Cross_Beta_pred_model_ExtraTree_1.3.1.pickle'
        ]
    
    for file_name in files_list:
        file_path = os.path.join(current_dir,file_name)
        if not os.path.isfile(file_path):
            print(f"File {file_path}" +RED+ "\tis missing..."+RESET)
            is_pb = True
        else:
            print(f"File {file_path}"+GREEN+ "\tOK"+RESET)
    
    if is_pb:
        print("Some folder and/or file are missing, try redownload Cross-Beta RF",
              "predictor to fixe this problem...")
    
    # Check library installation
    dependencies = [
        {'dep_name':'numpy',
         'import_name':'numpy',
         'version':'1.24.2'},
        {'dep_name':'pandas',
         'import_name':'pandas',
         'version':'2.0.0'},
        {'dep_name':'matplotlib',
         'import_name':'matplotlib',
         'version':'3.7.1'},
        {'dep_name':'scipy',
         'import_name':'scipy',
         'version':'1.8.0'},
        {'dep_name':'scikit-learn',
         'import_name':'sklearn',
         'version':'1.3.1'},
    ]
    
    for dependency in dependencies:
        test_lib = check_library(dependency)
        if not test_lib:
            is_pb = True
    if test_lib == 'wrong version':
        print('Some dependencies are not installed in the correct version. Cross-Beta',
              'predictor may work but the results can be invalid')

    return is_pb
    
