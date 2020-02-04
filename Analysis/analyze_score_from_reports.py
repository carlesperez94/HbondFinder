import os
import glob
import pandas as pd

RESULTS_PATH = "/home/carles/project/Almirall/ITK/Crystals/results/PLUS2/*"
SEP = "\t"
KEY_HB = "438H"
HB_COLUMNS = ["438H", "436O", "438O", "438H-ps", "436O-ps", "438O-ps"]
SCORE_COLUMN = " Binding Energy"


def find_reports_in_folder(folder, report_pref="report_*"):
    reports = glob.glob(os.path.join(folder, "*", report_pref))
    return reports


def fill_dataframe_from_reports(reports_list, separator="\t"):
    """
    From a list of reports (from the same simulation) it creates a single panda's dataframe which contains all
    information. WARNING: it ignores the first row (initial pose)
    :param reports_list: list of report files.
    :param separator: separator character between columns
    :return: pandas dataframe with report's content
    """
    df = pd.DataFrame()
    for report in reports_list:
        df_new = pd.read_csv(report, sep=separator, header=0)
        df_new = df_new[1:]  # Ignore the first row to avoid repeated snapshots
        df = df.append(df_new)
    return df


def filter_if_not_value(dataframe, column, value):
    return dataframe[dataframe[column] >= value]


def compute_pseudo_hb(df_row, column_pseudo_name, pseudo_id="-ps"):
    """
    It computes the quantity of pseudo hydrogen bonds. To do so you must provide a row of a dataframe. Given a pseudo-id,
    which is a substring that identifies the counts that also contain pseudo hydrogen bonds from the others, it computes
    the difference between the normal hydrogen bonds column and the one with additional pseudo HB.
    :param df_row: row of dataframe
    :param column_pseudo_name: column name of the row that contains pseudo HB
    :param pseudo_id: substring to identify the pseudo HB columns in the row
    :return: counting of pseudo HB in the column
    """
    column_normal_hb = column_pseudo_name.replace(pseudo_id, "")
    count_hb_normal = df_row[column_normal_hb]
    count_hb_pseudo = df_row[column_pseudo_name]
    number_of_pseudo = count_hb_pseudo - count_hb_normal
    return number_of_pseudo


def write_warnings_of_dataframe_size(dataframe, n_cutoff, folder, hbond_cutoff, file_object):
    if len(dataframe.index) < n_cutoff:
        msg = "WARNING: {} size is lower than {} ({}) for {} HB!\n".format(folder, n_cutoff, len(dataframe.index),
                                                                           hbond_cutoff)
        print(msg)
        file_object.write(msg)


def compute_hb_score(dataframe, hb_columns, pseudo_id="-ps", new_col_name="hb_score"):
    dataframe = dataframe.reset_index()  # We need unique indexes
    for n, itrow in enumerate(dataframe.iterrows()):
        idx, row = itrow
        hb_score = 0
        for col in hb_columns:
            if pseudo_id not in col:
                hb_score += row[col]
            else:
                pseudo_counter = compute_pseudo_hb(row, col, pseudo_id)
                hb_score += pseudo_counter/2  # Each pseudo H bond counts half
        dataframe.loc[idx, new_col_name] = hb_score
    return dataframe


def main(results_path_pattern, hb_key, hb_columns, score_column=" Binding Energy", sep="\t", hb_cutoff=3.0,
         n_cutoff=25):
    df_out = pd.DataFrame(columns=["System", "{}_mean".format(score_column), "{}_Q1_mean".format(score_column),
                                   "HBcutoff"])
    results_folders = sorted(glob.glob(results_path_pattern))
    logf = open("errors.log", "w")  # Open a file to store errors
    if not results_folders:
        raise FileNotFoundError("Any report in the following path: '{}'. Check if the path is correct."
                                .format(results_path_pattern))
    for folder in results_folders:
        print(folder)
        reports = sorted(find_reports_in_folder(folder))
        df = fill_dataframe_from_reports(reports, separator=sep)
        df_filtered = filter_if_not_value(df, column=hb_key, value=1)
        try:
            df_with_score = compute_hb_score(df_filtered, hb_columns)
            df_good_poses = filter_if_not_value(df_with_score, column="hb_score", value=hb_cutoff)
            hb_new_cutoff = hb_cutoff
            while df_good_poses.empty:
                hb_new_cutoff -= 0.5
                df_good_poses = filter_if_not_value(df_with_score, column="hb_score", value=hb_new_cutoff)
            write_warnings_of_dataframe_size(df_good_poses, n_cutoff, folder, hb_new_cutoff, logf)
            score = df_good_poses[score_column].mean()
            df_q1 = df_good_poses[df_good_poses[score_column] < df_good_poses[score_column].quantile(0.25)]
            score_q1 = df_q1[score_column].mean()
            df = pd.DataFrame({"System": [folder], "{}_mean".format(score_column): [score],
                               "{}_Q1_mean".format(score_column): [score_q1], "HBcutoff": [hb_new_cutoff]})
            df_out = df_out.append(df)
        except KeyError:
            print("{} DOES NOT HAVE ANY KEY INTERACTION!".format(folder))
    logf.close()  # Closing the error file
    df_out.to_csv("result_hb_analysis.tsv", sep=sep, header=True, index=False)  # Writing a tsv file with all the results


main(RESULTS_PATH, KEY_HB, HB_COLUMNS, SCORE_COLUMN, SEP, hb_cutoff=3, n_cutoff=25)
