#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv 

class CreateCSV:
    '''
    Using the results dictionary from motif_locator.py, generate a CSV file containing the enzyme name, motif sequence, cut site, number of observations, and the starting position of all searched enzymes in the plasmid.
    '''
    def __init__(self, results: dict, output_file: str):
        self.results = results
        self.output_file = output_file
        self.headers = ['enzyme', 'motif', 'cut_site', 'observed_count', 'start_positions']

    def create_csv_output(self):
        '''
        Creates the CSV output file in the `results` directory.
        '''
        with open(self.output_file, "w") as file:
            writer = csv.writer(file)
            writer.writerow(self.headers)
            
            for enzyme, (motif, cut_site, count, positions) in self.results.items():
                writer.writerow([enzyme, motif, cut_site, count, positions])
