#!/usr/bin/bash
# Change uid and gid here
sshfs -o reconnect,uid=1002,gid=1002 ubuntu@aws.big-data-biology.org:/share/work/Celio/AMPSphere/v2022_03/ data/aws_original_data

rsync --recursive --exclude data/aws_original_data/AMPSphere_generation_v.2022-03/analysis/helical_wheels* data/aws_original_data/* data/original_data
