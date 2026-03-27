import matplotlib.pyplot as plt
from pprint import pprint
import scipy.io
import numpy as np
from pathlib import Path
import spikeinterface as si  # import core only
import spikeinterface.extractors as se
import spikeinterface.preprocessing as spre
import spikeinterface.sorters as ss
import spikeinterface.postprocessing as spost
import spikeinterface.qualitymetrics as sqm
import spikeinterface.comparison as sc
import spikeinterface.exporters as sexp
import spikeinterface.curation as scur
import spikeinterface.widgets as sw
from pathlib import Path
import re
import os
import probeinterface as pi
import probeinterface.plotting as pi_plot
import cfad 
import json

global_job_kwargs = dict(n_jobs=4, chunk_duration="1s")
si.set_global_job_kwargs(**global_job_kwargs)


data = {
    "Behavioural": {
        "Aquisition": {
            "Extracted Behaviours": {
                "Darting" : {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Acquisition\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Acquisition\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Acquisition\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Acquisition\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Acquisition\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Acquisition\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Acquisition\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Acquisition\Behavioural Data\Extracted Behaviours\Darting",
                },
                "Freezing": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Acquisition\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Acquisition\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Acquisition\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Acquisition\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Acquisition\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Acquisition\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Acquisition\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Acquisition\Behavioural Data\Extracted Behaviours\Freezing",
                },
                "Grooming": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Acquisition\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Acquisition\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Acquisition\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Acquisition\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Acquisition\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Acquisition\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Acquisition\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Acquisition\Behavioural Data\Extracted Behaviours\Grooming",
                },
                "Position": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Acquisition\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Acquisition\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Acquisition\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Acquisition\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Acquisition\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Acquisition\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Acquisition\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Acquisition\Behavioural Data\Extracted Behaviours\Position",
                },
                "Rearing": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Acquisition\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Acquisition\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Acquisition\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Acquisition\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Acquisition\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Acquisition\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Acquisition\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Acquisition\Behavioural Data\Extracted Behaviours\Rearing",
                },
                "Velocity": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Acquisition\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Acquisition\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Acquisition\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Acquisition\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Acquisition\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Acquisition\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Acquisition\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Acquisition\Behavioural Data\Extracted Behaviours\Velocity",
                },
            },
            "Raw DLC Data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Acquisition\Behavioural Data\Raw DLC Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Acquisition\Behavioural Data\Raw DLC Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Acquisition\Behavioural Data\Raw DLC Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Acquisition\Behavioural Data\Raw DLC Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Acquisition\Behavioural Data\Raw DLC Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Acquisition\Behavioural Data\Raw DLC Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Acquisition\Behavioural Data\Raw DLC Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Acquisition\Behavioural Data\Raw DLC Data",
            },
            "SMM fitted data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Acquisition\Behavioural Data\SSM Fitted Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Acquisition\Behavioural Data\SSM Fitted Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Acquisition\Behavioural Data\SSM Fitted Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Acquisition\Behavioural Data\SSM Fitted Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Acquisition\Behavioural Data\SSM Fitted Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Acquisition\Behavioural Data\SSM Fitted Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Acquisition\Behavioural Data\SSM Fitted Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Acquisition\Behavioural Data\SSM Fitted Data",
            },
            "Traingulated Data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Acquisition\Behavioural Data\Triangulated Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Acquisition\Behavioural Data\Triangulated Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Acquisition\Behavioural Data\Triangulated Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Acquisition\Behavioural Data\Triangulated Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Acquisition\Behavioural Data\Triangulated Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Acquisition\Behavioural Data\Triangulated Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Acquisition\Behavioural Data\Triangulated Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Acquisition\Behavioural Data\Triangulated Data",
            },
            "Video Data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Acquisition\Behavioural Data\Video Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Acquisition\Behavioural Data\Video Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Acquisition\Behavioural Data\Video Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Acquisition\Behavioural Data\Video Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Acquisition\Behavioural Data\Video Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Acquisition\Behavioural Data\Video Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Acquisition\Behavioural Data\Video Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Acquisition\Behavioural Data\Video Data",
            },
        },
        "Extinction": {
            "Extracted Behaviours": {
                "Darting" : {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Behavioural Data\Extracted Behaviours\Darting",
                },
                "Freezing": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Behavioural Data\Extracted Behaviours\Freezing",
                },
                "Grooming": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Behavioural Data\Extracted Behaviours\Grooming",
                },
                "Position": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Behavioural Data\Extracted Behaviours\Position",
                },
                "Rearing": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Behavioural Data\Extracted Behaviours\Rearing",
                },
                "Velocity": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Behavioural Data\Extracted Behaviours\Velocity",
                },
            },
            "Raw DLC Data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Behavioural Data\Raw DLC Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Behavioural Data\Raw DLC Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Behavioural Data\Raw DLC Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Behavioural Data\Raw DLC Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Behavioural Data\Raw DLC Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Behavioural Data\Raw DLC Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Behavioural Data\Raw DLC Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Behavioural Data\Raw DLC Data",
            },
            "SMM fitted data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Behavioural Data\SSM Fitted Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Behavioural Data\SSM Fitted Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Behavioural Data\SSM Fitted Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Behavioural Data\SSM Fitted Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Behavioural Data\SSM Fitted Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Behavioural Data\SSM Fitted Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Behavioural Data\SSM Fitted Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Behavioural Data\SSM Fitted Data",
            },
            "Traingulated Data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Behavioural Data\Triangulated Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Behavioural Data\Triangulated Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Behavioural Data\Triangulated Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Behavioural Data\Triangulated Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Behavioural Data\Triangulated Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Behavioural Data\Triangulated Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Behavioural Data\Triangulated Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Behavioural Data\Triangulated Data",
            },
            "Video Data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Behavioural Data\Video Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Behavioural Data\Video Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Behavioural Data\Video Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Behavioural Data\Video Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Behavioural Data\Video Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Behavioural Data\Video Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Behavioural Data\Video Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Behavioural Data\Video Data",
            },
        
        },
        "Renewal": {
            "Extracted Behaviours": {
                "Darting" : {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Behavioural Data\Extracted Behaviours\Darting",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Behavioural Data\Extracted Behaviours\Darting",
                },
                "Freezing": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Behavioural Data\Extracted Behaviours\Freezing",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Behavioural Data\Extracted Behaviours\Freezing",
                },
                "Grooming": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Behavioural Data\Extracted Behaviours\Grooming",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Behavioural Data\Extracted Behaviours\Grooming",
                },
                "Position": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Behavioural Data\Extracted Behaviours\Position",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Behavioural Data\Extracted Behaviours\Position",
                },
                "Rearing": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Behavioural Data\Extracted Behaviours\Rearing",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Behavioural Data\Extracted Behaviours\Rearing",
                },
                "Velocity": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Behavioural Data\Extracted Behaviours\Velocity",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Behavioural Data\Extracted Behaviours\Velocity",
                },
            },
            "Raw DLC Data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Behavioural Data\Raw DLC Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Behavioural Data\Raw DLC Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Behavioural Data\Raw DLC Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Behavioural Data\Raw DLC Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Behavioural Data\Raw DLC Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Behavioural Data\Raw DLC Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Behavioural Data\Raw DLC Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Behavioural Data\Raw DLC Data",
            },
            "SMM fitted data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Behavioural Data\SSM Fitted Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Behavioural Data\SSM Fitted Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Behavioural Data\SSM Fitted Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Behavioural Data\SSM Fitted Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Behavioural Data\SSM Fitted Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Behavioural Data\SSM Fitted Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Behavioural Data\SSM Fitted Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Behavioural Data\SSM Fitted Data",
            },
            "Traingulated Data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Behavioural Data\Triangulated Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Behavioural Data\Triangulated Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Behavioural Data\Triangulated Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Behavioural Data\Triangulated Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Behavioural Data\Triangulated Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Behavioural Data\Triangulated Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Behavioural Data\Triangulated Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Behavioural Data\Triangulated Data",
            },
            "Video Data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Behavioural Data\Video Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Behavioural Data\Video Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Behavioural Data\Video Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Behavioural Data\Video Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Behavioural Data\Video Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Behavioural Data\Video Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Behavioural Data\Video Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Behavioural Data\Video Data",
            },
        
        },
    },
    "Combined": {
        "Aquisition": {
            "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Acquisition\Combined Data",
            "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Acquisition\Combined Data",
            "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Acquisition\Combined Data",
            "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Acquisition\Combined Data",
            "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Acquisition\Combined Data",
            "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Acquisition\Combined Data",
            "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Acquisition\Combined Data",
            "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Acquisition\Combined Data",
        },
        "Extinction": {
            "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Combined Data",
            "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Combined Data",
            "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Combined Data",
            "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Combined Data",
            "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Combined Data",
            "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Combined Data",
            "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Combined Data",
            "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Combined Data",
        },
        "Renewal": {
            "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Combined Data",
            "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Combined Data",
            "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Combined Data",
            "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Combined Data",
            "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Combined Data",
            "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Combined Data",
            "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Combined Data",
            "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Combined Data",
        },
    },
    "Neural": {
        "Aquisition": {
            "Concatenated Data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Acquisition\Neural Data\Concatenated Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Acquisition\Neural Data\Concatenated Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Acquisition\Neural Data\Concatenated Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Acquisition\Neural Data\Concatenated Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Acquisition\Neural Data\Concatenated Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Acquisition\Neural Data\Concatenated Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Acquisition\Neural Data\Concatenated Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Acquisition\Neural Data\Concatenated Data",
                "kilosort4": {
                    "Mouse 2": r"",
                    "Mouse 3": r"",
                    "Mouse 4": r"",
                    "Mouse 5": r"",
                    "Mouse 6": r"",
                    "Mouse 7": r"",
                    "Mouse 8": r"",
                    "Mouse 9": r"",
                },
            },
            "Raw Data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Acquisition\Neural Data\Raw Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Acquisition\Neural Data\Raw Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Acquisition\Neural Data\Raw Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Acquisition\Neural Data\Raw Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Acquisition\Neural Data\Raw Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Acquisition\Neural Data\Raw Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Acquisition\Neural Data\Raw Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Acquisition\Neural Data\Raw Data",
            },
            "Triggers": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Acquisition\Neural Data\Triggers",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Acquisition\Neural Data\Triggers",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Acquisition\Neural Data\Triggers",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Acquisition\Neural Data\Triggers",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Acquisition\Neural Data\Triggers",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Acquisition\Neural Data\Triggers",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Acquisition\Neural Data\Triggers",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Acquisition\Neural Data\Triggers",
            },
        },
        "Extinction": {
            "Concatenated Data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Neural Data\Concatenated Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Neural Data\Concatenated Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Neural Data\Concatenated Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Neural Data\Concatenated Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Neural Data\Concatenated Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Neural Data\Concatenated Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Neural Data\Concatenated Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Neural Data\Concatenated Data",
                "kilosort4": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Neural Data\Concatenated Data\kilosort4",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Neural Data\Concatenated Data\kilosort4",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Neural Data\Concatenated Data\kilosort4",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Neural Data\Concatenated Data\kilosort4",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Neural Data\Concatenated Data\kilosort4",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Neural Data\Concatenated Data\kilosort4",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Neural Data\Concatenated Data\kilosort4",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Neural Data\Concatenated Data\kilosort4",
                },
            },
            "Raw Data": {
                "Habituation": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Neural Data\Raw Data\mouse2_extinction_habituation_2025-06-06_13-13-19",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Neural Data\Raw Data\mouse3_extinction_habituation",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Neural Data\Raw Data\mouse4_extinction_habituation",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Neural Data\Raw Data\mouse5_extinction_habituation",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Neural Data\Raw Data\mouse6_extinction_habituation",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Neural Data\Raw Data\mouse7_extinction_habituation",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Neural Data\Raw Data\mouse8_extinction_habituation",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Neural Data\Raw Data\mouse9_extinction_habituation",
                },
                "part 1": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Neural Data\Raw Data\mouse2_extinction_p1_2025-06-06_13-24-00",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Neural Data\Raw Data\mouse3_extinction_p1",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Neural Data\Raw Data\mouse4_extinction_p1",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Neural Data\Raw Data\mouse5_extinction_p1",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Neural Data\Raw Data\mouse6_extinction_p1",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Neural Data\Raw Data\mouse7_extinction_p1",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Neural Data\Raw Data\mouse8_extinction_p1",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Neural Data\Raw Data\mouse9_extinction_p1",
                },
                "part 2": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Neural Data\Raw Data\mouse2_extinction_p2_2025-06-06_13-48-21",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Neural Data\Raw Data\mouse3_extinction_p2",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Neural Data\Raw Data\mouse4_extinction_p2",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Neural Data\Raw Data\mouse5_extinction_p2",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Neural Data\Raw Data\mouse6_extinction_p2",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Neural Data\Raw Data\mouse7_extinction_p2",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Neural Data\Raw Data\mouse8_extinction_p2",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Neural Data\Raw Data\mouse9_extinction_p2",
                },
            },
            "Triggers": {
                "Extracted events": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Neural Data\Triggers\mouse 2_extinction_extracted_events.mat",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Neural Data\Triggers\mouse 3_extinction_extracted_events.mat",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Neural Data\Triggers\mouse 4_extinction_extracted_events.mat",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Neural Data\Triggers\mouse 5_extinction_extracted_events.mat",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Neural Data\Triggers\mouse 6_extinction_extracted_events.mat",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Neural Data\Triggers\mouse 7_extinction_extracted_events.mat",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Neural Data\Triggers\mouse 8_extinction_extracted_events.mat",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Neural Data\Triggers\mouse 9_extinction_extracted_events.mat",
                },
                "Habituation": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Neural Data\Triggers\mouse2_extinction_habituation_triggers.mat",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Neural Data\Triggers\mouse3_extinction_habituation_triggers.mat",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Neural Data\Triggers\mouse4_extinction_habituation_triggers.mat",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Neural Data\Triggers\mouse5_extinction_habituation_triggers.mat",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Neural Data\Triggers\mouse6_extinction_habituation_triggers.mat",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Neural Data\Triggers\mouse7_extinction_habituation_triggers.mat",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Neural Data\Triggers\mouse8_extinction_habituation_triggers.mat",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Neural Data\Triggers\mouse9_extinction_habituation_triggers.mat",
                },
                "part 1": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Neural Data\Triggers\mouse2_extinction_p1_triggers.mat",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Neural Data\Triggers\mouse3_extinction_p1_triggers.mat",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Neural Data\Triggers\mouse4_extinction_p1_triggers.mat",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Neural Data\Triggers\mouse5_extinction_p1_triggers.mat",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Neural Data\Triggers\mouse6_extinction_p1_triggers.mat",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Neural Data\Triggers\mouse7_extinction_p1_triggers.mat",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Neural Data\Triggers\mouse8_extinction_p1_triggers.mat",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Neural Data\Triggers\mouse9_extinction_p1_triggers.mat",
                },
                "part 2": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Extinction\Neural Data\Triggers\mouse2_extinction_p2_triggers.mat",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Neural Data\Triggers\mouse3_extinction_p2_triggers.mat",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Neural Data\Triggers\mouse4_extinction_p2_triggers.mat",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Extinction\Neural Data\Triggers\mouse5_extinction_p2_triggers.mat",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Extinction\Neural Data\Triggers\mouse6_extinction_p2_triggers.mat",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Extinction\Neural Data\Triggers\mouse7_extinction_p2_triggers.mat",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Extinction\Neural Data\Triggers\mouse8_extinction_p2_triggers.mat",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Extinction\Neural Data\Triggers\mouse9_extinction_p2_triggers.mat",
                },
            },
        },
        "Renewal": {
            "Concatenated Data": {
                "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Neural Data\Concatenated Data",
                "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Neural Data\Concatenated Data",
                "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Neural Data\Concatenated Data",
                "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Neural Data\Concatenated Data",
                "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Neural Data\Concatenated Data",
                "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Neural Data\Concatenated Data",
                "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Neural Data\Concatenated Data",
                "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Neural Data\Concatenated Data",
                "kilosort4": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Neural Data\Concatenated Data\kilosort4",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Neural Data\Concatenated Data\kilosort4",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Neural Data\Concatenated Data\kilosort4",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Neural Data\Concatenated Data\kilosort4",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Neural Data\Concatenated Data\kilosort4",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Neural Data\Concatenated Data\kilosort4",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Neural Data\Concatenated Data\kilosort4",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Neural Data\Concatenated Data\kilosort4",
                },
            },
            "Raw Data": {
                "Checkerboard": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Neural Data\Raw Data\mouse2_renewal_checkerboard_2025-06-12_12-52-28",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Neural Data\Raw Data\mouse3_renewal_checkerboard",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Neural Data\Raw Data\mouse4_renewal_checkerboard",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Neural Data\Raw Data\mouse5_renewal_checkerboard",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Neural Data\Raw Data\mouse6_renewal_checkerboard",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Neural Data\Raw Data\mouse7_renewal_checkerboard",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Neural Data\Raw Data\mouse8_renewal_checkerboard",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Neural Data\Raw Data\mouse9_renewal_checkerboard",
                },
                "Habituation": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Neural Data\Raw Data\mouse2_renewal_habituation_2025-06-12_11-56-40",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Neural Data\Raw Data\mouse3_renewal_habituation",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Neural Data\Raw Data\mouse4_renewal_habituation",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Neural Data\Raw Data\mouse5_renewal_habituation",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Neural Data\Raw Data\mouse6_renewal_habituation",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Neural Data\Raw Data\mouse7_renewal_habituation",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Neural Data\Raw Data\mouse8_renewal_habituation",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Neural Data\Raw Data\mouse9_renewal_habituation",
                },
                "part 1": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Neural Data\Raw Data\mouse2_renewal_p1_2025-06-12_12-04-58",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Neural Data\Raw Data\mouse3_renewal_p1",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Neural Data\Raw Data\mouse4_renewal_p1",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Neural Data\Raw Data\mouse5_renewal_p1",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Neural Data\Raw Data\mouse6_renewal_p1",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Neural Data\Raw Data\mouse7_renewal_p1",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Neural Data\Raw Data\mouse8_renewal_p1",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Neural Data\Raw Data\mouse9_renewal_p1",
                },
                "part 2": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Neural Data\Raw Data\mouse2_renewal_p2_2025-06-12_12-28-14",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Neural Data\Raw Data\mouse3_renewal_p2",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Neural Data\Raw Data\mouse4_renewal_p2",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Neural Data\Raw Data\mouse5_renewal_p2",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Neural Data\Raw Data\mouse6_renewal_p2",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Neural Data\Raw Data\mouse7_renewal_p2",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Neural Data\Raw Data\mouse8_renewal_p2",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Neural Data\Raw Data\mouse9_renewal_p2",
                },
            },
            "Triggers": {
                "Checkerboard": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Neural Data\Triggers\mouse2_renewal_checkerboard_triggers.mat",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Neural Data\Triggers\mouse3_renewal_checkerboard_triggers.mat",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Neural Data\Triggers\mouse4_renewal_checkerboard_triggers.mat",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Neural Data\Triggers\mouse5_renewal_checkerboard_triggers.mat",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Neural Data\Triggers\mouse6_renewal_checkerboard_triggers.mat",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Neural Data\Triggers\mouse7_renewal_checkerboard_triggers.mat",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Neural Data\Triggers\mouse8_renewal_checkerboard_triggers.mat",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Neural Data\Triggers\mouse9_renewal_checkerboard_triggers.mat",
                },
                "Extracted events": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Neural Data\Triggers\mouse 2_renewal_extracted_events.mat",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Neural Data\Triggers\mouse 3_renewal_extracted_events.mat",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Neural Data\Triggers\mouse 4_renewal_extracted_events.mat",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Neural Data\Triggers\mouse 5_renewal_extracted_events.mat",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Neural Data\Triggers\mouse 6_renewal_extracted_events.mat",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Neural Data\Triggers\mouse 7_renewal_extracted_events.mat",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Neural Data\Triggers\mouse 8_renewal_extracted_events.mat",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Neural Data\Triggers\mouse 9_renewal_extracted_events.mat",
                },
                "Habituation": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Neural Data\Triggers\mouse2_renewal_habituation_triggers.mat",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Neural Data\Triggers\mouse3_renewal_habituation_triggers.mat",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Neural Data\Triggers\mouse4_renewal_habituation_triggers.mat",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Neural Data\Triggers\mouse5_renewal_habituation_triggers.mat",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Neural Data\Triggers\mouse6_renewal_habituation_triggers.mat",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Neural Data\Triggers\mouse7_renewal_habituation_triggers.mat",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Neural Data\Triggers\mouse8_renewal_habituation_triggers.mat",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Neural Data\Triggers\mouse9_renewal_habituation_triggers.mat",
                },
                "part 1": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Neural Data\Triggers\mouse2_renewal_p1_triggers.mat",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Neural Data\Triggers\mouse3_renewal_p1_triggers.mat",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Neural Data\Triggers\mouse4_renewal_p1_triggers.mat",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Neural Data\Triggers\mouse5_renewal_p1_triggers.mat",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Neural Data\Triggers\mouse6_renewal_p1_triggers.mat",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Neural Data\Triggers\mouse7_renewal_p1_triggers.mat",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Neural Data\Triggers\mouse8_renewal_p1_triggers.mat",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Neural Data\Triggers\mouse9_renewal_p1_triggers.mat",
                },
                "part 2": {
                    "Mouse 2": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 2\Renewal\Neural Data\Triggers\mouse2_renewal_p2_triggers.mat",
                    "Mouse 3": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Renewal\Neural Data\Triggers\mouse3_renewal_p2_triggers.mat",
                    "Mouse 4": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Renewal\Neural Data\Triggers\mouse4_renewal_p2_triggers.mat",
                    "Mouse 5": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 5\Renewal\Neural Data\Triggers\mouse5_renewal_p2_triggers.mat",
                    "Mouse 6": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 6\Renewal\Neural Data\Triggers\mouse6_renewal_p2_triggers.mat",
                    "Mouse 7": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 7\Renewal\Neural Data\Triggers\mouse7_renewal_p2_triggers.mat",
                    "Mouse 8": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Neural Data\Triggers\mouse8_renewal_p2_triggers.mat",
                    "Mouse 9": r"Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 9\Renewal\Neural Data\Triggers\mouse9_renewal_p2_triggers.mat",
                },
            },
        },
    }
}

def load_gain_from_oebin(recording_folder):
    """
    Parses structure.oebin to find the bit_volts (gain).
    Accepts the folder path and automatically appends 'structure.oebin'.
    """
    # 1. Construct the full file path
    oebin_path = os.path.join(recording_folder, "structure.oebin")
    
    # 2. Check if the FILE exists (not just the folder)
    if not os.path.isfile(oebin_path):
        # Fallback: check one level up (sometimes it sits in 'experiment1')
        parent_dir = os.path.dirname(recording_folder)
        parent_oebin = os.path.join(parent_dir, "structure.oebin")
        if os.path.isfile(parent_oebin):
            oebin_path = parent_oebin
        else:
            raise FileNotFoundError(f"structure.oebin not found at: {oebin_path}")

    # 3. Load and Parse
    with open(oebin_path, 'r') as f:
        meta = json.load(f)
    
    # Look for the Neural stream (usually has the most channels)
    # We loop through to find the one matching your channel count best, or just the largest
    neural_gain = None
    max_channels = 0
    
    for stream in meta.get('continuous', []):
        n_chans = stream.get('num_channels', 0)
        if n_chans > max_channels:
            max_channels = n_chans
            # Grab bit_volts from the first channel in the list
            if len(stream['channels']) > 0:
                neural_gain = stream['channels'][0]['bit_volts']
    
    if neural_gain is None:
        raise ValueError("Could not find a valid neural stream in structure.oebin")
    print(neural_gain)
    return neural_gain

def load_mice(start = 2, end = 9, rec = None, phase = None, slice = True):
    """
    Loads mice using a loop between start and end.
    If kilosorted then creates probe and assigns probe object from each mice
    Can build in to use raw ephys path to get probe info to avoid using ks
    returns lists for recordings, sortings and triggers
    fills initial elements of list with empty elements to line up iterations,
    i.e if loaded mice 2-3 then recordings will be [None, None, Mouse 2, Mouse 3]
    Done so the index needed to reference each mice is the same as mice number i.e recordings[2] = mice2 recording
    
    ToDo - 
    make so if no path provided then function call throws an exception for rec 
    smooth inputs for rec and phase to account for human error
    overall needs to smoothen the ui (data dict may need to be reworked)
    instead of array features need to orient to object oriented programming
    decouple trigger loading and slicing logic

    :param start: What mice to start at, default 2
    :param end: What mice to end at, default 9
    :param rec: type of recording (Habituation, part 1, part 2, or Checkerboard):
    :param phase: experimental phase, i.e what part of experiment it is rec from, Aqusition, Renewal or Extinction, can't find associated recording without it
    :param slice: Bool, default True, slice (and load) by triggers
    """
    recordings = []
    sortings = []
    triggers = []
    fs = 30000
    num_channels = 385
    dtype = "int16"
    buffer = 2 * fs
    


    for i in range(start):
        recordings.append(None)
        sortings.append(None)
        triggers.append(None)
    
    for i in range(start, end + 1):
        n = str(i)
        print("Mouse " + n)
        structure_path = data["Neural"][phase]["Raw Data"][rec]["Mouse " + n] + r"\Record Node 101\experiment1\recording1"
        gain = load_gain_from_oebin(structure_path)
        if slice:
            triggers_path = data["Neural"][phase]["Triggers"][rec]
            trigger = scipy.io.loadmat(triggers_path["Mouse " + n])
            trigger = np.array(trigger["evt"]).flatten().astype(int)
            #print(trigger)
            triggers.append(trigger)

        concat_kilosort = data["Neural"][phase]["Concatenated Data"]["kilosort4"]
        ks_path = concat_kilosort["Mouse " + str(i)]
        chan_pos_path = os.path.join(ks_path, "channel_positions.npy")
        chan_map_path = os.path.join(ks_path, "channel_map.npy")
        positions = np.load(chan_pos_path)
        chan_map = np.load(chan_map_path)
        chan_map = chan_map.flatten()

        probe = pi.Probe(ndim = 2, si_units ="um")
        probe.set_contacts(positions = positions, shapes = "square", shape_params = {"width": 12})
        probe.set_device_channel_indices(chan_map)
        probe.annotate(
            name='Neuropixels 2.0', 
            manufacturer='Imec', 
            description='Custom layout from ks folder'
        )


        rec_path = data["Neural"][phase]["Concatenated Data"]
        recording_path = rec_path["Mouse " + n] + r"\mouse" + n + r"_" + phase.lower() + r"_concatenated_neural_data.dat"
        recording = se.read_binary(file_paths = recording_path, dtype = dtype, sampling_frequency= fs, num_channels= num_channels, gain_to_uV= gain, offset_to_uV= 0)
        recording = recording.select_channels(channel_ids = np.arange(0, 384))
        recording= recording.frame_slice(start_frame = trigger[0] - buffer, end_frame = trigger[-1] + buffer)
        recording = recording.set_probe(probe, in_place=True)
        #print(recording)
        recordings.append(recording)

        sorting = se.read_kilosort(concat_kilosort["Mouse " + n])
        sorting = sorting.frame_slice(start_frame = trigger[0], end_frame = trigger[-1])
        sorting.register_recording(recording= recording)
        #print(sorting)
        sortings.append(sorting)

        print("Loaded")
        print("Triggers from " + str(trigger[0]) + " to " + str(trigger[-1]))
        print(recording)
        print("recording start at " + str(recording.get_start_time()) + " seconds")
        print("recording end at " + str(recording.get_end_time()) + " seconds")
        init_trigger_time = trigger[0] / fs
        last_trigger_time = trigger[-1] / fs
        print("Trigger start at " + str(init_trigger_time) + " seconds")
        print("Trigger final at " + str(last_trigger_time) + " seconds")
        print(sorting)
        print(probe)
        print("\n")

    return recordings, sortings, triggers


def export_ibl(recording, sorting, mark, n):
    """
    creates sorting analzyer and exports to ibl for specific recording and sorting
    mark denotes path mark, i.e checkerboard or habituation
    
    :param recording: recording object, best to pass a preprocessed recording
    :param sorting: sorting object
    :param mark: checkerboard or habituation, used to create output paths
    :param n: mouse number
    """
    sorting_analyzer = si.create_sorting_analyzer(sorting=sorting, recording=recording)

    # we need to compute some required extensions
    sorting_analyzer.compute(['random_spikes', 'templates', 'spike_amplitudes', 'spike_locations', 'noise_levels', 'quality_metrics'])
    # note that spike_locations are optional, but recommended to compute accurate spike depths

    # optionally, we can pass an LFP recording to compute RMS/PSD in the LFP band
    recording_lfp = spre.bandpass_filter(recording, freq_min=1, freq_max=300)
    # we can also decimate the LFP to speed up the process
    recording_lfp = spre.decimate(recording, 10)
    path = r"Z:\Saij\ephys\Mouse " + str(n) + mark

    sexp.export_to_ibl_gui(
    sorting_analyzer=sorting_analyzer,
    output_folder= path,
    lfp_recording=recording_lfp,
    n_jobs=-1
    )

    print("done")
