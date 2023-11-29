"""
Tips:
1. Do not rename the MyApp class and its five methods!
2. Return values of load_data, local_training, and global_aggregation methods will be communicated!
3. Return values of write_results will be discarded!
4. App supports five scenarios:
    * Native:
        - Centralized
        - Simulation of federated: No SMPC support
    * Containerized:
        - Centralized
        - Simulation of federated: No SMPC support
        - Federated
5. Containerized federated settings will be used in workflows for real world- collaboration
6. APP_NAME is case-sensitive and should contain the name of the app as a string and match the name of the config file
7. The config file should always named as `config.yml`
8. file structure:
    * Native
        - Centralized:
            data_dir to a directory containing all files for centralized training
        - Simulation:
9. transition loop between local_training and global_aggregation will continue until `last_round` attribute becomes True
"""

import utils

APP_NAME = "fc_example_config"


class MyAppFour(utils.AppFour):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.last_round = False

    def initial(self):
        return self.last_round

    def local_training(self, global_parameters):
        self.last_round = global_parameters
        return [None]

    def global_aggregation(self, local_parameters):
        self.last_round = True
        return self.last_round

    def write_results(self):
        pass

    def centralized(self):
        pass

