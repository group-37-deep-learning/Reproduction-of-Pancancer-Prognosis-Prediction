import pandas as pd
import json
import re
import numpy as np
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
# read the csv file
whole_excel = pd.read_csv('./data/pancancer_biospecimen.csv', sep= '\t')
# drop unnecessary columns
whole_excel = whole_excel.drop(whole_excel.columns[[0, 1, 2]], axis=1)
# drop unnecessary data from the cells
whole_excel['barcode'] = whole_excel['barcode'].map(lambda x: x[0:12])
whole_excel['project'] = whole_excel['project'].map(lambda x: x.replace("TCGA-", ""))
no_dup = whole_excel.drop_duplicates(subset=['barcode', 'project'])  # remove duplicates from the excel file
# print(no_dup)

# filter json data to the relevant vlaues
df = pd.read_json('./data/clinical.cases_selection.2020-03-08.json')
for i in range(0,len(df.demographic)):
    # get the necessary elements from the dataframe
    temp = [x for x in (str(df.demographic[i]).split(',')) if ("vital_status" in x) or
                                ("days_to_death" in x) or ("submitter_id" in x)]
    # result is a list, filter out unnecessary characters
    temp = list(filter(lambda x:x, map(lambda x:re.sub(r'[^A-Za-z-:123456789]', '', x), temp)))
    # only leave the TCGA_code and if alive
    # remove submitter id
    temp = [s.replace('submitterid:', '') for s in temp]
    # remove demographic
    temp = [s.replace('demographic', '') for s in temp]
    # remove vitalstatus
    temp = [s.replace('vitalstatus:', '') for s in temp]

    df.loc[i, 'demographic'] = temp


for i in range(1, len(df.diagnoses)):
    temp = [x for x in (str(df.diagnoses[i]).split(',')) if ("days_to_last_follow_up" in x)]
    temp = [s.replace('days_to_last_follow_up:', '')for s in temp]
    temp = list(filter(lambda x:x, map(lambda x:re.sub(r'[^1-9]', '', x), temp)))
    temp = ''.join(temp)
    temp = temp.replace(' ', '')
    df.loc[i, 'days_to_last_follow_up'] = temp

# delete unnecessary columns
del df['exposures']
del df['diagnoses']

# add new column to the dataframe
df['barcode'] = [''] * len(df.demographic)
df['days_to_death'] = [''] * len(df.demographic)
df['health'] = [''] * len(df.demographic)
# df['days_to_last_follow_up'] = [''] * len(df.demographic)
for i in range(0, len(df.demographic)):
    temp = str(df.demographic[i]).split(',') # split the element into parts
    for ii in range(0, len(temp)):
        if temp[ii].find('TCGA-')== 2:
            temp[ii] = temp[ii].replace("[", "")
            temp[ii] = temp[ii].replace("]", "")
            temp[ii] = temp[ii].replace("'", "")
            temp[ii] = temp[ii].replace(" ", "")
            df.loc[i,'barcode'] = temp[ii]
        elif temp[ii].find('daystodeath:')== 2:
            temp[ii] = temp[ii].replace("[", "")
            temp[ii] = temp[ii].replace("]", "")
            temp[ii] = temp[ii].replace("'", "")
            temp[ii] = temp[ii][12:len(temp[ii])]
            df.loc[i,'days_to_death'] = temp[ii]
        elif temp[ii].find('Alive')== 2:
            temp[ii] = temp[ii].replace("[", "")
            temp[ii] = temp[ii].replace("]", "")
            temp[ii] = temp[ii].replace("'", "")
            df.loc[i, 'health'] = temp[ii]
        elif temp[ii].find('Dead')== 2:
            temp[ii] = temp[ii].replace("[", "")
            temp[ii] = temp[ii].replace("]", "")
            temp[ii] = temp[ii].replace("'", "")
            df.loc[i, 'health'] = temp[ii]

jointed = pd.merge(df, no_dup, on='barcode', how='inner')
# there are no duplicates
# some have days_to_death = None remove those rows


jointed = jointed[jointed.days_to_death != 'None']


jointed =jointed.reset_index(drop=True) # reset index
# drop rows with faulty data, where the patient is alive but there is no days to last follow up
jointed.days_to_last_follow_up = jointed.days_to_last_follow_up.astype(str)
for i in range(0, len(jointed.health)):
    if (jointed.health[i] == ' Alive') and ((len(jointed.days_to_last_follow_up[i]) == 0 ) or
                                            (jointed.days_to_last_follow_up[i]== 'nan')):
        jointed = jointed.drop(i)
jointed = jointed.reset_index(drop=True)  # reset index
# create the plots

# cancer types
cancer_types = list(set(list(jointed.project)))
mortality = []
for i in range(0, len(cancer_types)):
    temp = jointed.loc[jointed['project'] == cancer_types[i]]
    total = len(temp)
    dead = len(temp.loc[jointed['health'] == ' Dead'])
    mortality.append(dead/ total)

# add the days_to_last_follow_up values to the days_to_death
for i in range(0, len(jointed.days_to_death)):
    if len(jointed.days_to_death[i]) == 0:
        jointed.days_to_death[i] = jointed.days_to_last_follow_up[i]
# check again remove empty rows
for i in range(0, len(jointed.days_to_death)):
    if len(jointed.days_to_death[i]) == 0:
        jointed = jointed.drop(i)
jointed = jointed.reset_index(drop=True)  # reset index

# sort the cancer types based on mortality values
mort_rate = [cancer_types for _,cancer_types in sorted(zip(mortality,cancer_types))]
#create the plots
# manually create the the list of the plots to look as similar as possible to the one in the article


left =['ACC', 'BRCA', 'CESC', 'DLBC', 'KICH', 'KIRC', 'KIRP', 'OV', 'PCPG', 'PRAD', 'SARC', 'SKCM', 'TGCT',
       'THCA', 'THYM', 'UCEC']
colors = ['crimson', 'green', 'gold', 'blue', 'orange', 'purple', 'lightblue', 'orchid', 'olive',
               'salmon', 'teal', 'pink', 'maroon', 'navy', 'silver', 'black']
right = ['BLCA', 'CHOL', 'COAD', 'ESCA', 'GBM', 'HNSC', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'PAAD', 'READ',
         'STAD', 'UCS', 'UVM']
kmf = KaplanMeierFitter()
ax = plt.subplot(1, 1, 1)
for i in range(0, len(left)):
    temp = jointed.loc[jointed['project'] == left[i]]
    days = temp.days_to_death.to_list()
    dead_or_not = []
    temp = temp.reset_index(drop=True)
    for ii in range(0, len(temp)):
        if temp.health[ii] == ' Alive':
            dead_or_not.append(0)
        else:
            dead_or_not.append(1)
    # convert the str lists to in strings
    days = list(map(int,days))
    kmf.fit(days, dead_or_not, label=left[i])
    ## Create an estimate
    ax = kmf.plot(ax=ax, ci_show=False, color=colors[i], linewidth=0.5)
    plt.legend(loc='upper center', ncol = 4, bbox_to_anchor=(.5, 1.5))
ax.set(xlim=(-500, 12500), ylim=(-0.1, 1.1))
# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.xlabel('Time')
plt.ylabel('Survival probability')
plt.savefig('left_plot21.png', dpi=1200,  bbox_inches='tight')
plt.show()


kmf = KaplanMeierFitter()
ax = plt.subplot(1, 1, 1)

for i in range(0, len(right)):
    temp = jointed.loc[jointed['project'] == right[i]]
    days = temp.days_to_death.to_list()
    dead_or_not = []
    #maxim = max(days)
    temp = temp.reset_index(drop=True)
    for ii in range(0, len(temp)):
        if temp.health[ii] == ' Alive':
            dead_or_not.append(0)
        else:
            dead_or_not.append(1)
    # convert the str lists to in strings
    days = list(map(int,days))
    kmf.fit(days, dead_or_not, label=right[i])
    ## Create an estimate
    ax = kmf.plot(ax=ax, ci_show=False, color=colors[i], linewidth=0.5)
    plt.legend(loc='upper center', ncol = 4, bbox_to_anchor=(.5, 1.5))
ax.set(xlim=(-500, 8500), ylim=(-0.1, 1.1))
# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.xlabel('Time')
plt.ylabel('Survival probability')
plt.savefig('right_plot21.png', dpi=1200,  bbox_inches='tight')
plt.show()

