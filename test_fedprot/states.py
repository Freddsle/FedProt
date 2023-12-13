from __future__ import annotations
from FeatureCloud.app.engine.app import AppState, app_state, Role, app
import time
import bios
import pandas as pd

CONFIG_FILE = "config.yml"
MOUNT_DIR = "/app"


@app_state(name='initial', role=Role.BOTH)
class InitialState(AppState):
    def register(self):
        self.register_transition('computing', label='Broadcast initial parameters')

    def run(self):
        self.store('aggregated', False)

        # read config
        self.read_config()
        self.log(f"Intensities: {self.load('intensities')}, Counts: {self.load('counts')}, Design: {self.load('design')}, Sep: {self.load('sep')}, Label: {self.load('label')}, Max NA Rate: {self.load('max_na_rate')}, Result Table: {self.load('result_table')}")

        # read data
        self.read_data()

        # initialize app
        self.log("Transition to computation...")
        return 'computing'

    def read_config(self):
        config = bios.read(f"{MOUNT_DIR}/{CONFIG_FILE}")
        if self.is_coordinator:
            self.store('client_name', "lab_A")
        else:
            self.store('client_name', "lab_B")
            
        self.store('intensities_name', config['intensities'])
        self.store('counts_name', config['counts'])
        self.store('design_name', config['design'])
        self.store('sep', config['sep'])
        self.store('label', config['label'])
        self.store('max_na_rate', config['max_na_rate'])
        self.store('result_table', config['result_table'])

    def read_data(self):
        self.log("Reading data...")
        intensities = pd.read_csv(f"{MOUNT_DIR}/data/{self.load('client_name')}/{self.load('intensities_name')}", sep=self.load('sep'), index_col=0)
        self.store('intensities_df', intensities)
        self.log(f"Data read! {intensities.shape[0]} samples, {intensities.shape[1]} proteins")


    def 



@app_state(name='computing')
class ComputingState(AppState):
    def register(self):
        self.register_transition('computing', role=Role.PARTICIPANT,
                                 label='Wait for another round of local training')
        self.register_transition('global_aggregation', role=Role.COORDINATOR,
                                 label='Collect and aggregate local models')
        self.register_transition('write_results', role=Role.BOTH, label='Broadcasting results')

    def run(self):
        self.log("Start computation...")
        time.sleep(3)
        self.log("Transition to aggregation or writing...")

        self.log(f"Aggregated: {self.load('aggregated')}")
        self.log(f"Coordinator: {self.is_coordinator}")
        
        if self.load('aggregated') or not self.is_coordinator:
            return 'write_results'
        else:
            return 'global_aggregation'
            

@app_state(name='global_aggregation')
class GlobalAggregation(AppState):
    def register(self):
        self.register_transition('computing', label='Broadcas final results')

    def run(self):
        self.log("Aggregating intermediate results...")
        self.store('aggregated', True)
        return 'computing'


@app_state(name='write_results')
class WriteResults(AppState):
    def register(self):
        self.register_transition('terminal', label='Done!')

    def run(self):
        self.log("Writing results...")
        self.log("Done!")
        return 'terminal'
