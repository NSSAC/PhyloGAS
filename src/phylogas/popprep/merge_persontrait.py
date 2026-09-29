import pandas as pd
import sys
import numpy as np

def calculate_smh_race(row):
    # default smh_race
    smh_race = 'Other'

    if int(row["hispanic_code"]) > 1:
        # Latino
        smh_race = 'Latino'
    elif int(row["race_code"]) == 1:
        # White
        smh_race = 'White'
    elif int(row["race_code"]) == 2:
        # Black
        smh_race = 'Black'
    elif int(row["race_code"]) == 6:
        # Asian
        smh_race = 'Asian'

    return smh_race

_, persontrait_file, person_file, fips_lookup_file, merged_file = sys.argv

pt_df = pd.read_csv(persontrait_file)
p_df = pd.read_csv(person_file)
fips_df = pd.read_csv(fips_lookup_file)

# Don't want to deal with duplicate columns
# 1.9.0 persontrait has pid,hid,age,age_group,gender,county_fips,home_latitude,home_longitude,admin1,admin2,admin3,admin4
# 1.9.0 person file has hid,pid,serialno,person_number,record_type,age,relationship,sex,school_enrollment,grade_level_attending,employment_status,occupation_socp,race,hispanic,designation
# Andrew wants: sex, age, home location for now, gender, race, hispanic, smh_race

p_df = p_df[[ "pid", "race","hispanic" ]]
p_df.set_index("pid")
pt_df.set_index("pid")

p_df.rename(columns={"race": "race_code", "hispanic": "hispanic_code"}, inplace=True)

pt_df.rename(columns={'age_group': 'age_group_code'}, inplace=True)

# merge the two dataframes
pt_df = pt_df.merge(p_df, on="pid", how="inner")

# add smh_race column
pt_df["smh_race"] = pt_df.apply(calculate_smh_race,axis=1)

# add the hispanic column
pt_df["latino"] = np.where(pt_df["hispanic_code"] > 1, True, False)

# add the race value
race_lookup = {1: "White alone", 2: "Black alone", 3: "Native American alone", 4: "Alaska Native alone", 5: "Native American and Alaska Native, tribe specified", 6: "Asian alone", 7: "Native Hawaiian and Other Pacific Islander alone", 8: "Some other race alone", 9: "Two or more races"}

pt_df["race"] = pt_df["race_code"].map(race_lookup)

age_group_lookup = {"p": "Preschool (0-4)", "s": "Student (5-17)", "a": "Adult (18-49)", "o": "Older adult (50-64)", "g": "Senior (65+)"}
pt_df["age_group"] = pt_df["age_group_code"].map(age_group_lookup)

# merge with fips code lookup
pt_df = pt_df.merge(fips_df, left_on="county_fips", right_on="FIPS", how="inner")

# remove unneeded columns
pt_df = pt_df.drop(["race_code","hispanic_code","age_group_code", "county_fips", "FIPS"], axis=1)

# write output to merged_file
pt_df.to_csv(merged_file, index=False)

