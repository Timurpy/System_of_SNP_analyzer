import pandas as pd
import yaml

# depricated
def parse_from_config(config_path="config.yaml"):
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    return config


def process_df_columns_ang_get_patient_list(df):
    columns = list(df.columns)
    renamed_columns = []
    unique_patients = []
    for patient_name in columns:
        if '(' in patient_name and ')' in patient_name:
            patient_name = patient_name.replace('(', '_v').replace(')', '')
        if '_v' in patient_name:
            patient_name = patient_name.split('_v')
            patient_id = int(patient_name[0].split()[-1])
            patient_probe_number = patient_name[1].replace('.', '_')
            new_name = f"patient_{patient_id}_{patient_probe_number}"
            renamed_columns.append(new_name)
            unique_patients.append(f"patient_{patient_id}")
        else:
            patient_id = int(patient_name.split()[-1])
            patient_probe_number = 0
            new_name = f"patient_{patient_id}_{patient_probe_number}"
            renamed_columns.append(new_name)
            unique_patients.append(f"patient_{patient_id}")

    unique_patients = list(set(unique_patients))
    df.columns = renamed_columns
    return unique_patients

def get_aggregated_patient_data_from_excel(path_to_input_file,
                                           seq_threshhold_depth,
                                           aggregation_method,
                                           ):
    # config = parse_from_config(config_path)
    
    # seq_threshhold_depth = config["seq_threshhold_depth"]
    # path_to_input_file = config["path_to_input_file"]
    # aggregation_method = config["aggregation_method"]
    if not isinstance(seq_threshhold_depth, int):
        try:
            seq_threshhold_depth = int(seq_threshhold_depth)
        except ValueError:
            print('AHTUNG!')
            seq_threshhold_depth = 10

    to_eng = {'Среднее': 'mean',
              'Минимальное': 'min',
              'Максимальное': 'max',
    }
    aggregation_method = to_eng[aggregation_method]
    # loading input excel table
    df = pd.read_excel(path_to_input_file, sheet_name=0)
    # preprocessing dataframe, casting it into convinience format for further analysis 
    df.replace({
        'yes, ': '', 
        'no': 0
        }, 
                regex=True, 
                inplace=True
                )
    df.fillna(0, inplace=True)
    df.index = list(map(lambda x: x.split()[-1], df['SNP']))
    df.drop(['SNP'], axis=1, inplace=True)
    df = df.astype(int)

    unique_patients = process_df_columns_ang_get_patient_list(df)

    for patient in unique_patients:
        df_patient = df[[col for col in df.columns if f"{patient}_" in col]]
        if aggregation_method == 'mean':
            df[patient] = df_patient.mean(axis=1).astype(int)
        elif aggregation_method == 'min':
            df[patient] = df_patient.min(axis=1).astype(int)
        elif aggregation_method == 'max':
            df[patient] = df_patient.max(axis=1).astype(int)
    
    #thresholding
    df[df < seq_threshhold_depth] = 0
            
    return df[unique_patients], unique_patients