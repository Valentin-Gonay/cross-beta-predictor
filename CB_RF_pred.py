''' 
Valentin Gonay
20/09/23
'''


# EXTERN IMPORTS 
import pickle # load the sci-kit learn mahine learning model
import os.path as path
from os import listdir
import re
import time
import argparse
import multiprocessing
import numpy as np


# INTERN IMPORTS
from utils import utils as utl
from utils import fold_pred as fpred
from utils import get_features
from utils import progress_bar as prg
from utils import check_install as ci

# Make prediction using the random forest model #
# Prediction for 1 sequence
def make_pred_1_seq(
        model, 
        sequence: str, 
        feature_list: list, 
        classification: str, 
        fold_pred_dict: dict
        ):
    '''Make the prediction for a sequence using a defined model with the given features and the 
    given classification method
    
    :param model: The sci-kit learn extratrees model that will give the prediction
    :type model: ExtraTreesClassifier

    :param sequence: The protein sequence to classify
    :type sequence: str 

    :param feature_list: The List of feature for the classification
    :type feature_list: list 

    :param classification: The classification method for converting AA into a group value. 
    Can be: 'mode_1', 'mode_2' or 'mode_3'
    :type classification: str 

    :param fold_pred_dict: The fold prediction provide by IUPred for the specific region
    :type fold_pred_dict: dict 

    :return: The list of predicted label (0 for Soluble IDR or 1 for Amyloid region) and The list \
    of prediction score with in each index, the score for Soluble IDR [i][0] and the score for \
    Amyloid [i][1] (ex: [[0.38 0.62]])
    :rtype: (list, list)
    '''
    
    entry = [
                get_features.get_features(
                    sequence, 
                    feature_list, 
                    classification
                    ) + [fold_pred_dict['iupred_region_pred']] 
            ]
    y_pred = model.predict(entry)[0]
    confidence = model.predict_proba(entry)[0]
    return y_pred,confidence


def run_prediction_fasta_longseq(
        model, 
        sequence_list: list, 
        feature_list: list, 
        classification: str, 
        pred_threshold: float = 0.5, 
        length_threshold: int | None = None
        ):
    '''Apply the prediction on a large sequence by cutting it in small fragment. Give a score for 
    each fragment. The final score for each amino acid correspond to the mean confidence.
    
    :param model: The sci-kit learn extratrees model that will give the prediction
    :type model: ExtraTreesClassifier

    :param sequence_list: The list of protein sequence you want to predict and will be fragmented. 
    Each entry must be a dictionary with a 'label' and a 'sequence' keys.
    :type sequence_list: list 

    :param feature_list: The List of feature for the classification
    :type feature_list: list 

    :param classification: The classification method for converting AA into a group value. 
    Can be: 'mode_1', 'mode_2' or 'mode_3'
    :type classification: str 

    :param pred_threshold: The threshold for the classification, every value higher than this 
    threshold will be considered as positive. If not provided, defaults = 0.5
    :type pred_threshold: float 

    :param length_threshold: The window size for the prediction. If not provided, defaults to None
    :type length_threshold: int | None 

    :return: The list of result of all the sequence of the given list
    :rtype: list
    '''

    result_list = []

    prg.print_loading_bar.start_time = prg.time.time() # Init start time
    iter = 0
    for entry in sequence_list:
        if len(sequence_list) > 1:
            prg.print_loading_bar(iter, len(sequence_list)-1, prefix='Computing prediction')
        entry_name = entry['label']
        entry_seq = entry['sequence']

        # Pass if the sequence don't match the 20 essential amino acids
        if re.match(r"^[ARNDCQEGHILKMFPSTWYV]*$", entry_seq): 
            if length_threshold == None:
                length_threshold = get_features.get_length_threshold(entry_seq)
                    
            amino_acid_list = []

            # Get folding prediction on total protein sequence
            total_seq_fold_score_dict = fpred.get_IUPred_allprot(entry_seq)

            # Run all the sequence, create fragment and make prediction on them            
            # The number of processors available
            process_count = multiprocessing.cpu_count() 

            # spread sequence fragment for in each processors
            sequence_fragment = get_seq_fragments(entry_seq, length_threshold)
            chunks = split_fragment_process(sequence_fragment, process_count)

            # Create sets of argument for each chunk (fairly spread with each available processors)
            arguments_lists = []
            for i in range(len(chunks)):
                arguments_lists.append((
                                    chunks[i], 
                                    total_seq_fold_score_dict, 
                                    classification, 
                                    feature_list, 
                                    model))

            # Make the prediction for each chunk in paralelle
            with multiprocessing.Pool(processes=process_count) as pool:
                results = pool.starmap(pred_seqFragment, arguments_lists)
            
            # Merge all results
            merged_amino_acid_list = []
            for result in results:
                merged_amino_acid_list += result[0]
            amino_acid_list = merge_AA_results(merged_amino_acid_list)
            
            # get the mean value for all the amino acids
            mean_list = []
            for i in range(len(amino_acid_list)):
                mean = np.mean(amino_acid_list[i]['score_list'])
                amino_acid_list[i]['mean_confidence'] = mean
                # extract and store in a new list, all the confidence mean
                temp_dict = {}
                temp_dict[amino_acid_list[i]['amino_acid']] = mean
                mean_list.append(temp_dict)
            
            region_pred_score_dict = {
                'iupred_region_pred': np.mean(total_seq_fold_score_dict['iupred_score_list'])
                }

            one_seq_result = {
                'prot_name': entry_name,
                'All_sequence_pred': make_pred_1_seq(
                                            model, 
                                            entry_seq, 
                                            feature_list, 
                                            classification, 
                                            fold_pred_dict=region_pred_score_dict
                                            )[1][1],
                'AA_list': amino_acid_list,
                'mean_list': mean_list,
                'AR_list': get_features.find_amyloid_region(mean_list, threshold = pred_threshold)
            }
            result_list.append(one_seq_result)
        iter += 1
    return result_list



# Multiprocessing specific functions
def get_seq_fragments(
        sequence: str, 
        window_size: int
        ):
    '''Generate fragment of the input sequence based on the window size, store the result in an 
    array where each element is a dict containing the sequence and the index of the AA in it
    
    :param sequence: The input amino-acid sequence
    :type sequence: str 

    :param window_size: The size of the window and the size of the generated fragment
    :type window_size: int 

    :return:The list of all the generated fragment where each list element is a dictionary 
    containing 'sequence' (the fragment sequence) and 'index' (the fragment AA index)
    :rtype: list
    '''

    result = []
    for i in range(len(sequence)-(window_size-1)):
        window = [i,i+window_size]
        fragment_seq = sequence[window[0]:window[1]]
        fragment_index = list(range(window[0],window[1]))
        result.append({'sequence':fragment_seq,'index':fragment_index})
    return result



def split_fragment_process(
        fragments: list, 
        nb_process: int
        ):
    '''Split the given sequence fragments list into chunks. 
    The number of chunks depend of the given number of processors
    
    :param fragments: The list of element to split
    :type fragments: list 

    :param nb_process: The number of processors available for the job
    :type nb_process: int 

    :return: The list of chunks based on the given number of process, each element of the list 
    is a list of fragments
    :rtype: list
    '''

    chunk_size = len(fragments) // nb_process
    if chunk_size == 0:
        chunk_size = 1
    chunks = [fragments[i:i+chunk_size] for i in range(0, len(fragments), chunk_size)]
    return chunks



def pred_seqFragment(
        fragments: list, 
        fold_pred: dict, 
        classification: str, 
        feature_list: list, 
        model
        ):
    '''Make the prediction on the given fragments with the given parameters
    
    :param fragments: The list of sequence fragment to predict
    :type fragments: list 

    :param fold_pred: The score for folding prediciton (ESM and IUPRed)
    :type fold_pred: dict 

    :param classification: The classification method
    :type classification: str 

    :param feature_list: The list of features
    :type feature_list: list 
    
    :param model: The sci-kit learn extratrees model that will give the prediction
    :type model: ExtraTreesClassifier

    :return: The Array with in index 0, the prediction result and in index 1, 
    the execution time for each fragments
    :rtype: list
    '''
    
    amino_acid_list = []
    time_exec_list = []
    for fragment in fragments:
        start_1_frg = time.time()
        for j in range(len(fragment['sequence'])):
            # check if the aminoacid already exist in the list
            if not any(d['index'] == fragment['index'][j] for d in amino_acid_list): 
                one_aminoacid_dict = {
                    'index':fragment['index'][j],
                    'amino_acid':fragment['sequence'][j],
                    'score_list':[],
                    'mean_confidence':0
                }
                amino_acid_list.append(one_aminoacid_dict)
        region_pred_score_dict = {
            'iupred_region_pred': np.mean(
                        fold_pred['iupred_score_list'][fragment['index'][0]:fragment['index'][-1]]
                        ),
        }
        y_pred, confidence = make_pred_1_seq(
                                    model, 
                                    fragment['sequence'], 
                                    feature_list, 
                                    classification, 
                                    fold_pred_dict=region_pred_score_dict
                                    )
        amyloid_confidence = confidence[1]
        
        # complete amino acids dicts with confidence values
        for aa_index in fragment['index']:
            for y in range(len(amino_acid_list)):
                if amino_acid_list[y]['index'] == aa_index:
                    amino_acid_list[y]['score_list'].append(amyloid_confidence)
        end_1_frg = time.time()
        total_1_frg = end_1_frg - start_1_frg
        time_exec_list.append(total_1_frg)
    return [amino_acid_list, time_exec_list]



def find_aaToMerge(aa_list: list):
    '''Find the amino acid that need to merge results
    
    :param aa_list: The array where every amino acid is a dictionary with a 'index' key
    :type aa_list: list

    :return: The list of AA to merge
    :rtype: list
    '''

    result = []
    treated_index = []
    for aa in aa_list:
        index_aa = aa['index']
        merge_aa = []
        if index_aa not in treated_index:
            treated_index.append(index_aa)
            for aa2 in aa_list:
                index_aa2 = aa2['index']
                if index_aa == index_aa2:
                    merge_aa.append(aa2)
        if len(merge_aa) > 0:
            result.append(merge_aa)
    return result



def merge_AA_results(amino_acid_list: list):
    '''Take a list of result by amino acids and merge the result of amino acids sharing the 
    same index
    
    :param amino_acid_list: The array of dictionaries where every entry corresponds to a amino acid 
    with its scores and info
    :type amino_acid_list: list 

    :return: The merge result in the same form that the given array of dictionaries
    :rtype: list
    '''

    to_merge_list = find_aaToMerge(amino_acid_list)
    merged_results = []
    for toMerge_grp in to_merge_list:
        aa_final = None
        for aa in toMerge_grp:
            if aa_final == None:
                aa_final = aa
            else:
                for score in aa['score_list']:
                    aa_final['score_list'].append(score) 
        merged_results.append(aa_final)
    return merged_results



# Main functions
def Cross_Beta_RF_pred(
        source: str, 
        source_type: str, 
        classification_method: str = 'mode_3', 
        threshold: float = 0.5, 
        label_col: str | None = None, 
        sequence_col: str | None = None, 
        window_size: int | None = None
        ):
    ''' Make the amyloidogenicity prediction for the input source 
    (csv, fasta file or sequence in String format) using specific treshold. 
    Give the result as a CSV file and a graph
    
    :param source: The path to access the csv or the fasta file. 
    Or the sequence without ID in String format.
    :type source: str 

    :param source_type: The type of the given source, must be 'csv', 'fastafile' or 'string'
    :type source_type: str 

    :param classification_method: The classification method for converting AA into a group value. 
    Can be: 'mode_1', 'mode_2' or 'mode_3'
    :type classification_method: str 

    :param threshold: The threshold use for the prediction, will determine the minimum confidence 
    score to predict a positive value. If not provided, defaults to 0.5
    :type threshold: float

    :param label_col: For csv files, give the column name where to find the sequence label. 
    If not provided, defaults to None.
    :type label_col: str | None

    :param sequence_col: For csv files, give the column name were to find the sequence. 
    If not provided, defaults to None.
    :type sequence_col: str | None

    :param window_size: The window size for the prediction of AR. If not provided, defaults to None.
    :type window_size: int | None

    :return: The list of result of all the sequence of the given list. Every result in the list is 
    a dictionary with all details about the sequence and the prediction
    :rtype: list
    '''

    ## Check if parameters are correct
    assert source_type in ['csv', 'fasta', 'sequence'], "Invalid format_type. Use 'csv', 'fasta', \
                                            or 'sequence'."
    assert threshold <= 1 and threshold >= 0, "Invalid threshold. Threshold value must be between \
                                            0 and 1."
    assert classification_method in ['mode_1', 'mode_2', 'mode_3'], "Invalid classification method.\
                                            Use 'mode_1', 'mode_2' or 'mode_3'."
    if source_type == 'csv':
        assert label_col != None or sequence_col != None, "Error, if the source type is csv, you \
                                            must give the label column name (label_col) and the \
                                            sequence column name (sequence_col)"
    
    # Get current dir
    current_directory = path.dirname(__file__)

    # Get model path
    relative_model_path = "data/Cross_Beta_pred_model_ExtraTree_1.3.1.pickle"
    model_path = path.join(current_directory, relative_model_path)

    ## load a trained sklearn random forest model
    forest = pickle.load(open(model_path, "rb"))

    feature_list = [
        'N', 'D', 'C', 'Q', 'I', 'L', 'M', 'F', 'T', 'W', 'Y', 
        'A_to_A', 'A_to_C', 
        'B_to_A', 'B_to_B', 'B_to_D', 'B_to_G', 
        'C_to_A', 'C_to_B', 'C_to_C', 'C_to_P', 'C_to_G', 
        'D_to_B', 'D_to_C', 'D_to_G', 
        'P_to_B', 'P_to_C', 'P_to_D', 'P_to_P', 'P_to_G', 
        'G_to_A', 'G_to_B', 'G_to_G', 
        'grp_A', 'grp_B', 'grp_C', 'grp_G', 
        'IUPred_score'
        ]

    ## Sequence extraction method
    if source_type == 'csv':
        sequence_list = utl.get_seq_from_csv(source, label_col, sequence_col)
    elif source_type == 'fasta':
        sequence_list = utl.get_seq_from_FASTA(source)
    elif source_type == 'sequence':
        sequence_list = utl.get_seq_from_strfasta(source)
    
    # Multi sequence prediction    
    result_list = run_prediction_fasta_longseq(
        forest, 
        sequence_list, 
        feature_list, 
        classification_method, 
        pred_threshold = threshold, 
        length_threshold = window_size
        )
    return(result_list)




### -------------------------------------------------------------------- ###
### ------------------------------- MAIN ------------------------------- ###
### -------------------------------------------------------------------- ###

if __name__ == '__main__':

    # ---------------------------------------------------------------------- #
    # ---------------------------- Get all args ---------------------------- #

    ### Create ArgumentParser object
    parser = argparse.ArgumentParser(
        description='Parsed argument for the usage of the Cross-Beta RF pred in commande line'
        )

    ### Add arguments
    ## Check install (optional)
    parser.add_argument(
        '-ci', 
        '--check_install', 
        action="store_true", 
        help='Check if all folder, files and python library necessary to run the \
            Cross-beta RF pred. are installed'
        )

    ## Input args
    # Input file:
    parser.add_argument(
        '-i',
        '--input', 
        type=str, 
        help='Input file path, can be file or folder'
        )

    # Input type:
    parser.add_argument(
        '-it',
        '--input_type', 
        type=str, 
        default='fasta', 
        help="Define the input type (sequence, csv or fasta) (default: 'fasta')."
        )

    ## For csv:
    # Name/id column
    parser.add_argument(
        '-nc',
        '--name_col', 
        type=str, 
        help="The name of the column containing the sequence ids. Only needed if \
            the input type is csv."
        )
    # Sequence column
    parser.add_argument(
        '-sc',
        '--seq_col', 
        type=str, 
        help="The name of the column containing the sequences. Only needed if the \
            input type is csv."
        )

    ## Other parameters
    # Classification method -> model trained for a classification method: 'mode_3', 
    # other classification methods are inplemented but may not provide accurate results
    # parser.add_argument(
    #   '-cl', 
    #   '--classification', 
    #   type=str, 
    #   default='mode_3', 
    #   help="Classification method, give the method use to group amino acid together \
    #       based on their characteristics"
    # )

    # Threshold (optional):
    parser.add_argument(
        '-t',
        '--threshold', 
        type=float, 
        default=0.54, 
        help="Classification threshold. Must be between 0 and 1 (Default: 0.54)"
        )

    # Predict window size
    parser.add_argument(
        '-ws', 
        '--prediction_window_size', 
        type=int, 
        default=0, 
        help="Specify the size of the window for the prediction. By default (if = 0) \
            the window size will adapt from 15 to 50 depending on the sequence length \
            (default: 0)"
        )

    # Draw graphs
    parser.add_argument(
        '-g', 
        '--draw_graph', 
        action="store_true", 
        help="Active or not the creation of graph for each result in the given output path"
        )

    # output files (optional):
    parser.add_argument(
        '-o',
        '--output', 
        type=str, 
        default='prediction_result.csv', 
        help="Output file name, the output results will always be saved in the 'results/' \
            folder (default: prediction_result.csv)"
        )

   
    args = parser.parse_args()




    # --------------------------------------------------------------------- #
    # --------------------------- Check install --------------------------- #
    if args.check_install:
        install_pb = ci.check_install()
        if install_pb:
            print("Problem in installation, please update your libraries and make sure all files \
                  and folder are present in the right place and correctly named")
        else:
            print("No problem detected in the installation")
        exit()





    # --------------------------------------------------------------------- #
    # ----------------------- Check all args values ----------------------- #
    
    # input
    try:
        utl.check_args_None(args.input, "--input")
    except ValueError as e:
        print(f"Error: {str(e)}")

    # input type
    assert args.input_type in ['sequence', 'fasta', 'csv'], f"Error: {args.input_type} is invalid. \
        --input_type must be a valid value: 'sequence', 'fasta' or 'csv'."
    # for csv only
    if args.input_type == 'csv':
        # name/id/label column
        try:
            utl.check_args_None(args.name_col, "--name_col")
        except ValueError as e:
            print(f"Error: {str(e)}")

        # sequence column
        try:
            utl.check_args_None(args.seq_col, "--seq_col")
        except ValueError as e:
            print(f"Error: {str(e)}")
    
    
    # prediction threshold
    utl.check_args_None(args.threshold, "--threshold")
    if args.threshold < 0 or args.threshold > 1:
        ValueError(f"Error: --threshold must be a between 0 and 1.")

    # output name
    utl.check_args_None(args.output, "--output")
    if not '.csv' in args.output:
        args.output += '.csv'
    

    if args.prediction_window_size == 0:
        args.prediction_window_size = None



    # -------------------------------------------------------------------- #
    # ------------------------- Get running mode ------------------------- #

    if args.input_type != 'sequence':
        # Get running mode by looking at the input path
        run_mode = utl.get_running_mode(args.input)

        file_list = []
        if run_mode == None:
            exit()
        elif run_mode == 'file':
            file_list.append(args.input)
        else:
            for file in listdir(args.input):
                if '.csv' in file or '.fasta' in file:
                    file_list.append(path.join(args.input,file))






    # -------------------------------------------------------------------- #
    # -------------------------- Init variables -------------------------- #

    # Get current dir
    current_directory = path.dirname(__file__)
    
    # Get result path
    relative_result_path = "results/"
    result_path = path.join(current_directory, relative_result_path)
    csv_path = path.join(current_directory, relative_result_path + args.output)
       
    predict_res = []



    # -------------------------------------------------------------------- #
    # ---------------------------- Prediction ---------------------------- #
    
    start_time = time.time()

    print("\nProcessing data...")

    if args.input_type == 'sequence':
        predict_res = Cross_Beta_RF_pred(   args.input,
                                            source_type = 'sequence', 
                                            classification_method = 'mode_3', 
                                            threshold = args.threshold, 
                                            window_size = args.prediction_window_size
                                        )
        
        print("\nSaving results ...")
        # Display / save results
        utl.save_amino_acid_pred_result_csv(csv_path, predict_res)

        if args.draw_graph:
            for res in predict_res:
                utl.drow_graph_result(
                    res['mean_list'], 
                    plot_name = res['prot_name'], 
                    threshold_line = args.threshold, 
                    save = result_path
                    )

        print("Job done, Thank you for using Cross-Beta RF predictor.\n")



    elif args.input_type == 'fasta':
        for file in file_list:
            if '.fasta' in file:
                print(f"\nPrediction on file {file}...")
                result_list = Cross_Beta_RF_pred(   file,  
                                                    source_type = 'fasta',
                                                    classification_method = 'mode_3',
                                                    threshold = args.threshold, 
                                                    window_size = args.prediction_window_size 
                                                )
                for result in result_list:
                    predict_res.append(result)
        
        print("\nSaving results ...")
        # Display / save results
        utl.save_amino_acid_pred_result_csv(csv_path, predict_res)

        if args.draw_graph:
            for res in predict_res:
                utl.drow_graph_result(
                    res['mean_list'], 
                    plot_name = res['prot_name'], 
                    threshold_line = args.threshold, 
                    save = result_path
                    )

        print("Job done, Thank you for using Cross-Beta RF predictor.\n")


                
    else:
        for file in file_list:
            if '.csv' in file:
                print(f"\nPrediction on file {file}...")
                result_list = Cross_Beta_RF_pred(   file,
                                                    source_type = 'csv',
                                                    classification_method = 'mode_3',
                                                    threshold = args.threshold, 
                                                    label_col= args.name_col,  
                                                    sequence_col= args.seq_col, 
                                                    window_size = args.prediction_window_size 
                                                )

                for result in result_list:
                    predict_res.append(result)
        
        print("\nSaving results ...")
        # Display / save results
        utl.save_amino_acid_pred_result_csv(csv_path, predict_res)

        if args.draw_graph:
            for res in predict_res:
                utl.drow_graph_result(
                    res['mean_list'], 
                    plot_name = res['prot_name'], 
                    threshold_line = args.threshold, 
                    save = result_path
                    )

        print("Job done, Thank you for using Cross-Beta RF predictor.\n")





    end_time = time.time()
    duration_time = end_time - start_time
    
    # Format remaining time
    minutes = int(duration_time // 60)
    seconds = int(duration_time % 60)
    
    duration_time_formated = f'{minutes:02}:{seconds:02}'
    print(f'\nExecution time: {duration_time_formated}')

