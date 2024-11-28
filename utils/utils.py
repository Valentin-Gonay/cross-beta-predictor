''' 
Valentin Gonay
16/10/23
'''


import pandas as pd
import csv
import matplotlib.pyplot as plt
import os


### Utils function ###
# Extract sequence from set, DB or FASTA/Multi-FASTA file
def get_seq_from_setfile(file_path: str):
    '''Get all sequence and the corresponding label from a set file (must contains a 'LABEL', 
    a 'Cluster_ID' and a 'Sequence' columns)
    
    :param file_path: The path to access to the set file
    :type file_path: str 
    
    :return: An array of all the sequence and their label organized as a dictionnary \
    (ex: [\
        {'label': 'Amyloid', \
        'cluster':'42', \
        'sequence': 'ABEEGG'},\
        {'label': 'Soluble', \
        'cluster':'101', \
        'sequence':'GEDDAG'}\
        ]\
    )
    :rtype: list
    '''

    seq_array = []
    df = pd.read_csv(file_path,sep=';',header=0).fillna("NaN")
    for row in df.itertuples():
        entry_dict = {'label':'', 'cluster':row.Cluster_ID, 'sequence':row.Sequence}
        if row.LABEL == 'AMYLOID':
            entry_dict['label'] = 'Amyloid'
            seq_array.append(entry_dict)
        else:
            entry_dict['label'] = 'Soluble'
            seq_array.append(entry_dict)
    return seq_array



def get_seq_from_DBfile(
        file_path: str, 
        sequence_colname: str = 'AR Sequence'
        ):
    '''Get all sequence and the corresponding label from a DB file (must contains a 'LABEL', 
    a 'Cluster_ID' and a sequence columns)
    
    :param file_path: The path to access to the set file
    :type file_path: str 

    :param sequence_colname: The name of the column containing the sequence to predict. 
    If not provided, defaults to 'AR Sequence'
    :type sequence_colname: str 
    
    :return: An array of all the sequence and their label organized as a dictionnary \
    (ex: [\
        {'label': 'Amyloid', \
        'cluster':'42', \
        'sequence': 'ABEEGG'},\
        {'label': 'Soluble', \
        'cluster':'101', \
        'sequence':'GEDDAG'}\
        ]\
    )
    :rtype: list
    '''
    
    seq_array = []
    df = pd.read_csv(file_path,sep=';',header=0).fillna("NaN")
    for index, row in df.iterrows():
        entry_dict = {'label':'', 'cluster':row['Cluster_ID'], 'sequence':row[sequence_colname]}
        if row['LABEL'] == 'AMYLOID':
            entry_dict['label'] = 'Amyloid'
            seq_array.append(entry_dict)
        else:
            entry_dict['label'] = 'Soluble'
            seq_array.append(entry_dict)
    return seq_array


def get_seq_from_FASTA(file_path: str):
    '''Get all sequence and the corresponding name from a FASTA file 
    
    :param file_path: The path to access to the FASTA file
    :type file_path: str 

    :return: An array of all the sequence and their label organized as a dictionnary 
    (ex: [{'label': 'protein_1', 'sequence': 'ABEEGG'},{'label': 'protein_2', 'sequence':'GEDDAG'}])
    :rtype: list
    '''

    sequences_list = []
    current_sequence = ''
    current_entry = {}
    
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()

        # Check if the first line starts with '>'
        if first_line.startswith('>'):
            current_entry['label'] = first_line[1:]
        else:
            current_entry['label'] = 'sequence_query'
            current_sequence = first_line


        for line in file:
            line = line.strip()
            # print('line_test: ',line)

            # New sequence
            if line.startswith('>'):
                # store the current sequence in the current entry and reset its value
                if current_sequence is not None:
                    current_entry['sequence'] = current_sequence
                    current_sequence = None

                # store the previous entry  in the list if a new sequence start and reset its value
                if current_entry != {}:
                    sequences_list.append(current_entry)
                    current_entry = {}
                
                current_entry['label'] = line[1:]
                current_sequence = ''
            
            # if the line corresponds to a new sequence line (same entry)
            else:
                current_sequence += line

        # Store the last line as the sequence of the last entry    
        if current_sequence is not None:
            current_entry['sequence'] = current_sequence
            sequences_list.append(current_entry)
            
    return sequences_list



def get_seq_from_csv(
        csv_file_path: str, 
        label_colname: str, 
        sequence_colname: str
        ):
    '''Convert a csv file into a python dict. Use specific column names to extract information 
    from the csv file
    
    :param csv_file_path: The path and file name of the source CSV file.
    :type csv_file_path: str 

    :param label_colname: The name of the column where to get the fasta entry name/identifier.
    :type label_colname: str 

    :param sequence_colname: The name of the column where to get the protein sequence
    :type sequence_colname: str 

    :return: The list of entry of the created fasta file. Each entry is a dictionary composed of a 
    'label' and a 'sequence' keys.
    :rtype: list
    '''

    # open CSV and extract informations
    entry_list = []
    df = pd.read_csv(csv_file_path,sep=';',header=0).fillna("NaN")
    for index, row in df.iterrows():
        entry_list.append({'label':row[label_colname], 'sequence':row[sequence_colname]})
    
    return entry_list



def get_seq_from_strfasta(sequences: str):
    '''Extract sequences from a fasta / multifasta in string format
    
    :param sequences: The input sequences in fasta or multifasta in String format 
    (not the path to the file)
    :type sequences: str 

    :return: An array of all the sequence and their label organized as a dictionnary 
    (ex: [{'label': 'protein_1', 'sequence': 'ABEEGG'},{'label': 'protein_2', 'sequence':'GEDDAG'}])
    :rtype: list
    '''

    sequences_list = []
    current_sequence = None
    current_entry = {}

    file = sequences.split('\n')
    for i in range(len(file)):
        line = file[i].strip()
        if i == 0 and not line.startswith('>'):
            current_entry['label'] = 'sequence_query'
            current_sequence = ''
        
        if line.startswith('>'):
            # store the current sequence in the current entry and reset its value
            if current_sequence is not None:
                current_entry['sequence'] = current_sequence
                current_sequence = None

            # store the previous entry in the list if a new sequence start and reset its value
            if current_entry != {}:
                sequences_list.append(current_entry)
                current_entry = {}

            current_entry['label'] = line[1:]
            current_sequence = ''

        # if the line corresponds to a new sequence line (same entry)
        else:
            current_sequence += line
        
    # Store the last line as the sequence of the last entry    
    if current_sequence is not None:
        current_entry['sequence'] = current_sequence
        sequences_list.append(current_entry)

    return sequences_list

        


# Display the result in the console
def display_res_1_pred(
        entry: dict, 
        pred: int, 
        confidence: float
        ):
    '''Display result of one prediction associated with the cluster ID with the target label is 
    there is one, the predicted label and the confidence score.
    
    :param entry: The entry used for the prediction with information about the label, the cluster 
    id and the sequence
    :type entry: dict 

    :param pred: The predicted label (0 for Soluble or 1 for Amyloid)
    :type pred: int 

    :param confidence: The confidence score of the prediction (from 0.5 to 1.0)
    :type confidence: float 

    :return: The dictionary containing all the displayed information (cluster, label, prediction, 
    confidence score) for the given entry
    :rtype: dict
    '''

    if pred == 1:
        prediction = 'Amyloid'
    else:
        prediction = 'Soluble'
    print(
        'Cluster:',entry['cluster'], 
        '\tLabel:', entry['label'], 
        '\tPrediction:', prediction,
        '\tConfidence score:', confidence)
    
    return {
        'cluster':entry['cluster'], 
        'label':entry['label'], 
        'prediction':prediction, 
        'confidence':confidence
        }



# Save the prediction result in a CSV file
def save_res_csv(
        save_path: str, 
        entries: list, 
        col_names: list
        ):
    '''Save result in CSV format of one prediction associated with the cluster ID with the target 
    label is there is one, the predicted label and the confidence score.
    
    :param save_path: The path and name for the creation of the result file
    :type save_path: str 

    :param entries: The list of entry containing all information about the entry and its 
    prediction results
    :type entries: list 

    :param col_names: The list of column names for the CSV file. All given column names 
    have to be Strings
    :type col_names: list 
    '''

    if not '.csv' in save_path:
        save_path += '.csv'
    
    with open(save_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=';')
        
        # Write header
        csv_writer.writerow(col_names)
        
        # Write data rows
        for entry in entries:
            csv_writer.writerow(
                [entry['cluster'], 
                 entry['label'], 
                 entry['prediction'], 
                 entry['confidence']])

    print(f'{save_path} created successfully')


def save_amino_acid_pred_result_csv(
        save_path: str, 
        result_list: list
        ):
    '''Save the prediction result giving amyloidogenicity score by amino acids in a csv file
    
    :param save_path: The path of the file and its name
    :type save_path: str 

    :param result_list: The list of result where every item is a dictionary containing 'prot_name', 
    'AA_list', 'All_sequence_pred', 'mean_list', and'AR_list' keys.  
    :type result_list: list 
    '''

    if not '.csv' in save_path:
        save_path += '.csv'

    colnames_list = [
        'Query_name',
        'Sequence_length', 
        'Average_protein_prediction', 
        'AR_position', 
        'Amino_acids_score'
        ]
    
    with open(save_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=';')

        # Write header
        csv_writer.writerow(colnames_list)
        
        # Write data rows
        for entry in result_list:
            csv_writer.writerow(
                [entry['prot_name'], 
                 len(entry['mean_list']), 
                 entry['All_sequence_pred'], 
                 entry['AR_list'], 
                 entry['mean_list']])

    print(f'{save_path} created successfully')





# Display the prediction result by amino acids as a bar chart
def drow_graph_result(
        aa_mean_confidence_list: list, 
        plot_name: str = '', 
        threshold_line: float = 0.5, 
        display: bool = False, 
        save: str = ''
        ):
    '''Drow a graph showing the mean confidence score for all the amino acids of the given list.
    
    :param aa_mean_confidence_list: The list of dictionary where the key is the 1 letter amino acid 
    code and the value is the mean confidence score
    :type aa_mean_confidence_list: list 

    :param plot_name: Specify the name the plot should take. If not provided, defaults to ''.
    :type plot_name: str 

    :param threshold_line:  Give the value where to drow the threshold line. The value must be 
    between 0 and 1. If 0, don't drow the line. If not provided, defaults to 0.5
    :type threshold_line: float 

    :param display: If True, display one by one all the graphic results. 
    If not provided, defaults to False. 
    :type display: bool 

    :param save: Give the specific file path where to save the created plots. 
    if save = None, don't save the plots. If not provided, defaults to ''.
    :type save: str 
    '''
    # clear plot in memory
    plt.clf()

    # extract values and labels in 2 lists
    labels = range(len(aa_mean_confidence_list))
    values = [list(d.values())[0] for d in aa_mean_confidence_list]
    
    
    # create bar char
    plt.bar(labels,values)

    if threshold_line != 0:
        plt.axhline(y = threshold_line, color='r', linestyle='-')
    
    plt.title(plot_name)
    plt.ylim(0,1)
    plt.xlabel('Sequence index')
    plt.ylabel('Amyloidogenicity score')

    plt.xticks(range(0, len(labels), 10))

    if display:
        plt.show()
    
    if save != None:
        plt.savefig(save+plot_name+'_result_graph.png', format='png')




def get_running_mode(input_path: str):
    '''Select the runing mode depending on the input path: folder or file
    
    :param input_path: The path of the input (file or folder)
    :type input_path: str 
    
    :return: The running mode: on folder: will try all files of the folder (only on .fasta or .csv) 
    or on one given file.
    :rtype: str
    '''

    if os.path.exists(input_path):
        if os.path.isfile(input_path):
            return 'file'
        elif os.path.isdir(input_path):
            return 'folder'
        else:
            print(f"{input_path} exists but is neither a file nor a directory.")
            return None
    else:
        print(f"{input_path} does not exist.")
        return None



def check_args_None(
        value, 
        args_name: str
        ):
    '''Return error message if the given value is None
    
    :param value: The given value, can be any type
    :type value: any
    :param args_name: The name of the parameter, used to edit the error message.
    :type args_name: str 
    '''
    
    if value is None:
        error_message = f"{args_name} is None, please provide a valid value for this argument."
        raise ValueError(error_message)
