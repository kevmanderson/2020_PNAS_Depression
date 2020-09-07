#!/bin/python

import os
import pandas as pd
import subprocess
import feather
import numpy as np
import calendar
import copy
from datetime import datetime
import scipy as sp


# this script will format and create data fields in the UKB dataframe
#    1. Rename UKB fields/dataframe columns
#    2. Derive probable MDD status
#    3. Neuroticism scale summed score
#    4. 


def force_symlink(src, dest):
    try:
        os.symlink(src, dest)
    except:
        os.remove(dest)
        os.symlink(src, dest)


def recode_yes(val):
    if pd.isna(val):
        return np.nan
    elif val == 'Yes':
        return 1
    else:
        return 0


def derive_probable_mdd(df, visit):
    '''
    Take the pandas data frame as input and derive MDD phenotypes

    (1) MDD/Anxiety treatment seeking
    (2) Probable single/moderate/severe depression

    :param df:
    :param visit:
    :return:
    '''

    # if they've ever sought professional treatment (Yes=1, No=0)
    treat_seek = np.where((df['f.2090.' + visit + '.0'] == 'Yes') | (df['f.2100.' + visit + '.0'] == 'Yes'))[0]
    treat_arr = df['f.2090.' + visit + '.0'].copy(deep=True)
    treat_arr.loc[treat_seek] = 'Yes'

    # create new variable indicating whether they've ever sought treatment
    df['mdd_treatment_seek.' + visit + '.0'] = treat_arr

    # probable single MDD
    # 4598 - Ever depressed for a whole week
    # 4609 - Longest period of depression
    # 4620 - Number of depression episodes
    depr_single_idxs = np.where((df['f.4598.' + visit + '.0'] == 'Yes') & (df['f.4609.' + visit + '.0'] >= 2) & (
            df['f.4620.' + visit + '.0'] == 1) & (treat_arr == 'Yes'))[0]

    # 4631 - Ever unenthusiastic/disinterested for a whole week
    # 5375 - ongest period of unenthusiasm / disinterest
    # 5386 - Number of unenthusiastic/disinterested episodes
    anhed_single_idxs = np.where((df['f.4631.' + visit + '.0'] == 'Yes') & (df['f.5375.' + visit + '.0'] >= 2) & (
            df['f.5386.' + visit + '.0'] == 1) & (treat_arr == 'Yes'))[0]

    # idxs of single episode subjects
    single_idxs = np.unique(np.concatenate([depr_single_idxs, anhed_single_idxs]))

    # create array of zeros (maintaining nan values)
    prob_single_episode = pd.Series([recode_yes(x) for x in df['f.4598.' + visit + '.0']])
    prob_single_episode[~np.isnan(prob_single_episode)] = 0

    # 1= single episode subjects
    prob_single_episode.iloc[single_idxs] = 1

    df['prob_single_MDD.' + visit + '.0'] = prob_single_episode

    # probable moderate MDD
    mod_depr_single_idxs = np.where((df['f.4598.' + visit + '.0'] == 'Yes') & (df['f.4609.' + visit + '.0'] >= 2) & (
            df['f.4620.' + visit + '.0'] >= 2) & (df['f.2090.' + visit + '.0'] == 'Yes'))[0]
    mod_anhed_single_idxs = np.where((df['f.4631.' + visit + '.0'] == 'Yes') & (df['f.5375.' + visit + '.0'] >= 2) & (
            df['f.5386.' + visit + '.0'] >= 2) & (df['f.2090.' + visit + '.0'] == 'Yes'))[0]

    # idxs of moderate episode subjects
    mod_idxs = np.unique(np.concatenate([mod_depr_single_idxs, mod_anhed_single_idxs]))

    # probable severe MDD
    sev_depr_single_idxs = np.where((df['f.4598.' + visit + '.0'] == 'Yes') & (df['f.4609.' + visit + '.0'] >= 2) & (
            df['f.4620.' + visit + '.0'] >= 2) & (df['f.2100.' + visit + '.0'] == 'Yes'))[0]
    sev_anhed_single_idxs = np.where((df['f.4631.' + visit + '.0'] == 'Yes') & (df['f.5375.' + visit + '.0'] >= 2) & (
            df['f.5386.' + visit + '.0'] >= 2) & (df['f.2100.' + visit + '.0'] == 'Yes'))[0]

    # idxs of severe episode subjects
    sev_idxs = np.unique(np.concatenate([sev_depr_single_idxs, sev_anhed_single_idxs]))

    # get rid of moderate subjecst that also meet criteria for severe
    new_mod_idxs = [x for x in mod_idxs if x not in sev_idxs]

    prob_mod_episode = pd.Series([recode_yes(x) for x in df['f.4598.' + visit + '.0']])
    prob_mod_episode[~np.isnan(prob_mod_episode)] = 0
    prob_mod_episode.iloc[new_mod_idxs] = 1
    df['prob_moderate_MDD.' + visit + '.0'] = prob_mod_episode

    prob_sev_episode = pd.Series([recode_yes(x) for x in df['f.4598.' + visit + '.0']])
    prob_sev_episode[~np.isnan(prob_sev_episode)] = 0
    prob_sev_episode.iloc[sev_idxs] = 1
    df['prob_severe_MDD.' + visit + '.0'] = prob_sev_episode

    mdd_arr = pd.Series([recode_yes(x) for x in df['f.4598.' + visit + '.0']])
    mdd_arr[~np.isnan(mdd_arr)] = 0
    mdd_arr.iloc[single_idxs] = 1
    mdd_arr.iloc[new_mod_idxs] = 2
    mdd_arr.iloc[sev_idxs] = 3
    df['prob_MDD.' + visit + '.0'] = mdd_arr

    return df


def sum_neuroticism(df, visit):
    '''
    Simple summation of the 12 neuroticism items

    :param df:
    :param visit:
    :return:
    '''
    neuroticism_fields = ['f.' + str(num) + '.' + visit + '.0' for num in np.arange(1920, 2040, 10)]
    neurot_mat = pd.DataFrame()
    for f in neuroticism_fields:
        f_num = pd.Series([recode_yes(x) for x in df[f]])
        neurot_mat = pd.concat([neurot_mat, f_num], axis=1)
    df['neuroticism.' + visit + '.0'] = neurot_mat.sum(axis=1, skipna=False)
    return df


def derive_bipolar(df, visit):
    '''
    Probable Bipolar I/II

    :param df:
    :param visit:
    :return:
    '''
    # if they've ever sought professional treatment (Yes=1, No=0)
    bp_cardinal = np.where((df['f.4642.' + visit + '.0'] == 'Yes') | (df['f.4653.' + visit + '.0'] == 'Yes'))[0]
    bp_cardinal_arr = df['f.4642.' + visit + '.0'].copy(deep=True)
    bp_cardinal_arr.loc[bp_cardinal] = 'Yes'

    symptoms = ['I was more active than usual', 'I was more talkative than usual', 'I needed less sleep than usual',
                'I was more creative or had more ideas than usual']
    symptom_mat = pd.DataFrame()
    for col in range(0, 4):
        symptom_present = pd.Series([1 if val in symptoms else 0 for val in df['f.6156.' + visit + '.' + str(col)]])
        symptom_mat = pd.concat([symptom_mat, symptom_present], axis=1)

    symptom_counts = symptom_mat.sum(axis=1, skipna=False)
    symptom_counts[np.where(df['f.6156.' + visit + '.0'] == 'All of the above')[0]] = 4

    bp1_single_idxs = np.where(
        (bp_cardinal_arr == 'Yes') & (symptom_counts >= 3) & (df['f.5663.' + visit + '.0'] == 'A week or more') & (df[
                                                                                                                       'f.5674.' + visit + '.0'] == 'Needed treatment or caused problems with work, relationships, finances, the law or other aspects of life'))[
        0]
    bp2_single_idxs = np.where(
        (bp_cardinal_arr == 'Yes') & (symptom_counts >= 3) & (df['f.5663.' + visit + '.0'] == 'A week or more') & (df[
                                                                                                                       'f.5674.' + visit + '.0'] != 'Needed treatment or caused problems with work, relationships, finances, the law or other aspects of life'))[
        0]

    bp_len_arr = np.zeros(len(df['f.5663.' + visit + '.0']))
    bp_len_arr[np.where(df['f.5663.' + visit + '.0'].isnull())[0]] = np.nan
    bp_len_arr[bp1_single_idxs] = 1
    bp_len_arr[bp2_single_idxs] = 2
    bp_len_arr = pd.Series(bp_len_arr)

    df['bipolar_status.' + visit + '.0'] = bp_len_arr
    return df


def combine_bp_mdd(df, visit):
    '''
    Take the inferred MDD/BP values and combine into a single variable

    :param df:
    :param visit:
    :return:
    '''
    combined_arr = df['prob_MDD.' + visit + '.0'].copy()
    combined_arr.loc[df['bipolar_status.' + visit + '.0'] == 1] = 4
    combined_arr.loc[df['bipolar_status.' + visit + '.0'] == 2] = 5
    df['bp_mdd_status.' + visit + '.0'] = combined_arr
    return df


def sum_depression(df, visit):
    '''
    Sum 7 binary depression items and 1 ordinal weight symptom

    :param df:
    :param visit:
    :return:
    '''
    sad = pd.Series([recode_yes(x) for x in df['f.20446.' + visit + '.0']])
    interest = pd.Series([recode_yes(x) for x in df['f.20441.' + visit + '.0']])
    sleep = pd.Series([recode_yes(x) for x in df['f.20532.' + visit + '.0']])
    worthless = pd.Series([recode_yes(x) for x in df['f.20450.' + visit + '.0']])
    conc = pd.Series([recode_yes(x) for x in df['f.20435.' + visit + '.0']])
    tired = pd.Series([recode_yes(x) for x in df['f.20449.' + visit + '.0']])
    suicide = pd.Series([recode_yes(x) for x in df['f.20437.' + visit + '.0']])

    weight = np.zeros(len(df['f.20536.' + visit + '.0']))
    weight[np.where(df['f.20536.' + visit + '.0'].isnull())[0]] = np.nan
    weight_idxs = np.where(
        (df['f.20536.' + visit + '.0'] == 'Lost weight') | (df['f.20536.' + visit + '.0'] == 'Gained weight') | (
                df['f.20536.' + visit + '.0'] == 'Both gained and lost some weight during the episode'))
    weight[weight_idxs] = 1
    weight = pd.Series(weight)

    count_df = pd.DataFrame(np.vstack([sad, interest, sleep, worthless, conc, tired, suicide, weight]))
    symptom_ct = count_df.sum()
    na_counts = count_df.isna().sum()
    symptom_ct[na_counts == 8] = np.nan

    df['mdd_online_symptom_sum.' + visit + '.0'] = symptom_ct
    return df


def recode_fields(df, field_df):
    '''
    Replace the UKB field ids with meaningful column names

    :param df:
    :param field_df:
    :return:
    '''
    new_col_arr = []
    for col in df.columns:
        try:
            splits = col.split('.')
            new_name = field_df.loc[field_df.id == int(splits[1]), 'name'].values[0]
            splits[1] = new_name
            new_col = '.'.join(splits[1:])
            print(col + ' >> ' + new_col)
        except:
            new_col = col
            print(col + ' >> ' + new_col)
        new_col_arr.append(new_col)
    df.columns = new_col_arr
    return df


def days_between(d1, d2):
    d1 = datetime.strptime(d1, "%Y-%m-%d")
    d2 = datetime.strptime(d2, "%Y-%m-%d")
    return abs((d2 - d1).days) / 365


def infer_age(df_recode):
    '''
    The age field at time of scan has lots of missing values. So calculate age at scan from birthday and scanning visit date.

    :param df:
    :return:
    '''
    yob = df_recode['year_of_birth.0.0'].astype(object).astype(str).str.replace('[.]0', '')
    mob = df_recode['month_of_birth.0.0']

    mo_dict = dict((v, k) for k, v in enumerate(calendar.month_name))
    mo_num = pd.Series([str(mo_dict[val]).zfill(2) if not pd.isnull(val) else np.nan for val in mob])
    yob_mo = yob.str.cat(mo_num, sep='-')
    yob_mo_day = yob_mo + '-01'
    scan_date = df_recode['Date_of_attending_assessment_centre.2.0']

    df_recode['age_at_scan'] = [
        days_between(str(yob_mo_day[k]), str(scan_date[k])) if not pd.isnull(scan_date[k]) else np.nan for i, k in
        enumerate(range(len(scan_date)))]
    df_recode['age_square'] = np.square(df_recode['age_at_scan'])
    return df_recode


def scz_bp_online_symptoms(df):
    # distress related to psychosis
    df['OMH_Psychosis_distress_ordinal'] = copy.deepcopy(df['f.20462.0.0'])
    mapping = {'Not distressing at all, it was a positive experience': 1,
               'Not distressing, a neutral experience': 2,
               'A bit distressing': 3,
               'Quite distressing': 4,
               'Very distressing': 5,
               'Do not know': np.nan,
               'Prefer not to answer': np.nan}
    df['OMH_Psychosis_distress_ordinal'] = df['OMH_Psychosis_distress_ordinal'].replace(mapping)

    # frequency of psychotic experiences
    df['OMH_Psychosis_freqExpr_lastYr_ordinal'] = copy.deepcopy(df['f.20467.0.0'])
    mapping = {'Not at all': 1,
               'Once or twice': 2,
               'Less than once a month': 3,
               'More than once a month': 4,
               'Nearly every day or daily': 5,
               'Prefer not to answer': np.nan}
    df['OMH_Psychosis_freqExpr_lastYr_ordinal'] = df['OMH_Psychosis_freqExpr_lastYr_ordinal'].replace(mapping)

    # length of manic episode
    df['OMH_Mania_lenManiaIrrit_ordinal'] = copy.deepcopy(df['f.20492.0.0'])
    mapping = {'Less than 24 hours': 1,
               'At least a day, but less than a week': 2,
               'A week or more': 3,
               'Prefer not to answer': np.nan,
               'Do not know': np.nan}
    df['OMH_Mania_lenManiaIrrit_ordinal'] = df['OMH_Mania_lenManiaIrrit_ordinal'].replace(mapping)

    # type of BP symptoms
    bp_strings = [x for x in pd.unique(df['f.20548.0.1'])[1:]]
    bp_strings = pd.Series(bp_strings).loc[pd.notna(bp_strings)].tolist()
    bp_cols = df.columns[df.columns.str.contains('20548')]
    for symptom in bp_strings:
        symptom_short = symptom.replace('I was more ', '').replace(' than usual', '')
        field_name = 'OMH_bipolar_' + symptom_short.replace(' ', '_') + '_0'
        df[field_name] = np.nan
        df.loc[pd.notnull(df['f.20548.0.1']), field_name] = 0
        for cur_col in bp_cols:
            df.loc[df[cur_col].astype(object) == symptom, field_name] = 1

    # symptom severity
    df['OMH_Mania_severity_num'] = copy.deepcopy(df['f.20493.0.0'])
    mapping = {'No problems': 0,
               'Needed treatment or caused problems with work, relationships, finances, the law or other aspects of life.': 1,
               'Prefer not to answer': np.nan,
               'Do not know': np.nan}
    df['OMH_Mania_severity_num'] = df['OMH_Mania_severity_num'].replace(mapping)
    return (df)


def mental_health_icd(df):
    icd_cols = df.columns[df.columns.str.contains('41204|41202')]

    # schizophrenia codes
    scz_codes = ['F20' + str(i) for i in np.arange(0, 10)]
    scz_codes = scz_codes + ['F210', 'F220', 'F228', 'F229', 'F230', 'F231', 'F233', 'F239', 'F250', 'F251', 'F252',
                             'F258', 'F259', 'F280', 'F290']
    df['icd_scz'] = np.nan
    df.loc[pd.notnull(df['f.41202.0.0']), 'icd_scz'] = 0
    for col in icd_cols:
        df.loc[df[col].isin(scz_codes), 'icd_scz'] = 1

    # bipolar codes
    bp_codes = ['F31' + str(i) for i in np.arange(0, 10)]
    bp_codes = bp_codes + ['F300', 'F301', 'F302', 'F308', 'F309']
    df['icd_bipolar'] = np.nan
    df.loc[pd.notnull(df['f.41202.0.0']), 'icd_bipolar'] = 0
    for col in icd_cols:
        df.loc[df[col].isin(bp_codes), 'icd_bipolar'] = 1

    # depression
    depr_codes = ['F320', 'F321', 'F322', 'F323', 'F328', 'F329', 'F330', 'F331', 'F332', 'F333', 'F334', 'F338',
                  'F339', 'F340', 'F341', 'F348', 'F380', 'F388', 'F390']
    df['icd_depression'] = np.nan
    df.loc[pd.notnull(df['f.41202.0.0']), 'icd_depression'] = 0
    for col in icd_cols:
        df.loc[df[col].isin(depr_codes), 'icd_depression'] = 1

    return (df)


def process_pheno(df_path, field_df_path):
    base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_sst'
    field_df_path = os.path.join(base_dir, 'reference_files/ukb_field_dict.csv')
    df_path = os.path.join(base_dir, 'data/ukb/enc/ukb_pheno_mri.csv')
    df_feather = os.path.join(base_dir, 'data/ukb/enc/ukb_pheno_mri.feather')

    field_df = pd.read_csv(field_df_path)
    print('> reading UKB dataframe')
    # df = pd.read_csv(df_path)
    df = feather.read_dataframe(df_feather)

    print('> infer single/moderate/severe MDD')
    df = derive_probable_mdd(df=df, visit='0')
    df = derive_probable_mdd(df=df, visit='1')
    df = derive_probable_mdd(df=df, visit='2')

    print('> infer bipolar I/II')
    df = derive_bipolar(df=df, visit='0')
    df = derive_bipolar(df=df, visit='1')
    df = derive_bipolar(df=df, visit='2')

    print('> create combined MDD/BP variable')
    df = combine_bp_mdd(df=df, visit='0')
    df = combine_bp_mdd(df=df, visit='1')
    df = combine_bp_mdd(df=df, visit='2')

    print('> recode online SCZ/BP symptoms')
    df = scz_bp_online_symptoms(df=df)

    print('> scz/bipolar/depr ICD diagnosis')
    df = mental_health_icd(df=df)

    print('> sum online depression inventory')
    df = sum_depression(df, visit='0')

    print('> calculate neuroticism scores')
    df = sum_neuroticism(df, visit='0')
    df = sum_neuroticism(df, visit='1')
    df = sum_neuroticism(df, visit='2')

    print('> calculate age at imaging visit')
    df_recode = recode_fields(df, field_df)
    df_recode = df_recode.loc[df_recode['year_of_birth.0.0'].notna()]
    df_recode = df_recode.reindex(np.arange(df_recode.shape[0]))
    df_recode = infer_age(df_recode)

    for visit in range(3):
        print('calc diastology/systolic BP:' + str(visit))
        df_recode['diastolic.' + str(visit) + '.0'] = pd.concat(
            [df_recode['diastolic_auto.' + str(visit) + '.0'], df_recode['diastolic_manual.' + str(visit) + '.0']],
            axis=1).mean(1)
        df_recode['systolic.' + str(visit) + '.0'] = pd.concat(
            [df_recode['systolic_auto.' + str(visit) + '.0'], df_recode['systolic_manual.' + str(visit) + '.0']],
            axis=1).mean(1)

    # convert treatment seeking to numeric
    df_recode['mdd_treatment_seek_num.2.0'] = df_recode['mdd_treatment_seek.2.0']
    mapping = {'Yes': 1, 'No': 0, 'Do not know': 0, 'Prefer not to answer': 0}
    df_recode = df_recode.replace({'mdd_treatment_seek_num.2.0': mapping})

    df_recode['mdd_status.2.0'] = copy.deepcopy(df_recode['bp_mdd_status.2.0'])
    df_recode.loc[df_recode['mdd_status.2.0'] == 4, 'mdd_status.2.0'] = np.nan
    df_recode.loc[df_recode['mdd_status.2.0'] == 5, 'mdd_status.2.0'] = np.nan

    for x in df_recode.columns[df_recode.columns.str.contains('\.x')]:
        print(x.replace('.x', ''))
        df_recode[x.replace('.x', '')] = df_recode[x]

    print('> remove subjects with mismatched genetic/reported sex')
    male_1 = np.where(df_recode['genetic_sex.0.0'] == 'Male')[0]
    male_2 = np.where(df_recode['sex.0.0'] == 'Male')[0]
    male_mismatch = np.setdiff1d(male_1, male_2)

    female_1 = np.where(df_recode['genetic_sex.0.0'] == 'Female')[0]
    female_2 = np.where(df_recode['sex.0.0'] == 'Female')[0]
    female_mismatch = np.setdiff1d(female_1, female_2)
    genetic_mismatch = np.sort(np.concatenate([male_mismatch, female_mismatch]))

    df_recode = df_recode.drop(genetic_mismatch)

    # get rid of some unnecessary columns
    df_recode = df_recode.drop(df_recode.columns[df_recode.columns.str.contains('cancer')], axis=1)
    df_recode = df_recode.drop(df_recode.columns[df_recode.columns.str.contains('illness_noncancer')], axis=1)
    df_recode = df_recode.drop(df_recode.columns[df_recode.columns.str.contains('illness_noncancer')], axis=1)
    df_recode = df_recode.drop(df_recode.columns[df_recode.columns.str.contains('\.y')], axis=1)
    df_recode = df_recode.drop(df_recode.columns[df_recode.columns.str.contains('\.x')], axis=1)
    df_recode = df_recode.drop(df_recode.columns[df_recode.columns.str.contains('lenFirstPress')], axis=1)
    df_recode = df_recode.drop(df_recode.columns[df_recode.columns.str.contains('numPress')], axis=1)
    df_recode = df_recode.drop(df_recode.columns[df_recode.columns.str.contains('dMRI')], axis=1)

    print('>> save processed feather-df')
    out_path = df_path.replace('ukb_pheno', 'ukb_pheno_processed')
    print(out_path)
    df_recode.to_csv(out_path, index=False)
    feather.write_dataframe(df_recode, df_path.replace('ukb_pheno.csv', 'ukb_pheno_processed.feather'))

    return (df_recode)
