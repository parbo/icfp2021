#!/bin/bash

for i in {1..132}
do
    curl -H "Authorization: Bearer $ICFP_API_KEY" -o problems/$i.problem https://poses.live/api/problems/$i
done
