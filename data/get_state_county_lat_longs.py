import pandas as pd
import sys

_, state_code, output_file = sys.argv


df = pd.read_csv("state_county_centroids.csv")

# State,FIPS,County,Latitude,Longitude

# Filter it down to state_code
df = df.loc[df["State"] == state_code]

# Convert to format
df["geo_resolution"] = "county"

df["county"] = df["County"] + "," + df["State"]

df = df [[ "geo_resolution","county","Latitude","Longitude" ]]

# write out, no index, no headers, to tsv
df.to_csv(output_file, index=False, header=False, sep='\t')

