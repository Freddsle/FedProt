"""
    FeatureCloud Four States App Template
    Copyright 2023 Mohammad Bakhtiari. All Rights Reserved.
    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at
        http://www.apache.org/licenses/LICENSE-2.0
    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
"""
from __future__ import annotations
from FeatureCloud.app.engine.app import AppState, app_state, Role, app


def create_client(app_four_instance, featurecloud_app_instance, centralized=False):
    """ Generate states for a client based on centralized or federated scenario

    Parameters
    ----------
    app_four_instance: MyApp
    featurecloud_app_instance: app
    centralized: bool

    Returns
    -------
    general_app_instance
    """
    if centralized:
        return generate_centralized_states(app_four_instance, featurecloud_app_instance)
    return generate_federated_states(app_four_instance, featurecloud_app_instance)


def generate_federated_states(app_four: "AppFour", featurecloud_app: None or app = None):
    '''The function `generate_federated_states` generates a set of states for a federated learning
    application using the FeatureCloud framework.
    
    Parameters
    ----------
    app_four : AppFour
        The parameter `app_four` is an instance of the `AppFour` class. It is used to access the methods
    and attributes of the `AppFour` class within the `generate_federated_states` function.
    featurecloud_app : None or app
        The `featurecloud_app` parameter is an instance of the FeatureCloud application. It is used to
    define the states and transitions of the federated learning process.
    
    Returns
    -------
        the `featurecloud_app` object.
    
    '''

    @app_state(name='initial', role=Role.BOTH, app_instance=featurecloud_app)
    class InitialState(AppState):
        def register(self):
            self.register_transition('Local_Training', label='Broadcast initial parameters')

        def run(self):
            data_to_broadcast = app_four.initial()
            if self.is_coordinator:
                self.broadcast_data(data_to_broadcast)
            elif data_to_broadcast is not None:
                RuntimeWarning("Participants data broadcast request was ignored!")

            return 'Local_Training'

    @app_state('Local_Training', app_instance=featurecloud_app)
    class LocalTraining(AppState):

        def register(self):
            self.register_transition('Global_Aggregation', role=Role.COORDINATOR,
                                     label='Collect and aggregate local models')
            self.register_transition('Local_Training', role=Role.PARTICIPANT,
                                     label='Wait for another round of local training')
            self.register_transition('Write_Results', role=Role.BOTH, label='Scape the loop')

        def run(self):
            received_data = self.await_data()
            data_to_send = app_four.local_training(global_parameters=received_data)
            if data_to_send:
                self.send_data_to_coordinator(data_to_send,
                                              use_smpc=app_four.config["use_smpc"])
            if app_four.last_round or not self.is_coordinator:
                return "Write_Results"
            if self.is_coordinator:
                return "Global_Aggregation"

    @app_state('Global_Aggregation', app_instance=featurecloud_app)
    class GlobalAggregation(AppState):

        def register(self):
            self.register_transition('Local_Training', role=Role.COORDINATOR,
                                     label='Broadcat the global models and go to the next next round')

        def run(self):
            if app_four.config['use_smpc']:
                local_data = self.aggregate_data()
            else:
                local_data = self.gather_data()

            data_to_broadcast = app_four.global_aggregation(local_parameters=local_data)

            self.log("Data to broadcast...")
            print(data_to_broadcast)

            self.broadcast_data(data_to_broadcast)
            return 'Local_Training'

    @app_state('Write_Results', app_instance=featurecloud_app)
    class WriteResults(AppState):

        def register(self):
            self.register_transition('terminal', label='Finish app execution')

        def run(self):
            app_four.write_results()
            return 'terminal'

    return featurecloud_app


def generate_centralized_states(app_four: "AppFour", featurecloud_app: app):
    '''The function generates centralized states for an app using the AppFour and app classes.
    
    Parameters
    ----------
    app_four : AppFour
        The parameter "app_four" is an instance of the class "AppFour". It is being used in the function
    "generate_centralized_states" to perform some operations related to centralized training.
    featurecloud_app : app
        The parameter `featurecloud_app` is an instance of the `app` class. It is used to define the states
    and transitions of the application.
    
    Returns
    -------
        the modified featurecloud_app object.
    
    '''

    @app_state(name='initial', role=Role.BOTH, app_instance=featurecloud_app)
    class InitialState(AppState):

        def register(self):
            self.register_transition('Centralized', label='Immediate transition')

        def run(self):
            return 'Centralized'

    @app_state('Centralized', app_instance=featurecloud_app)
    class Centralized(AppState):

        def register(self):
            self.register_transition('terminal', label='Finish centralized training')

        def run(self):
            app_four.centralized()
            return 'terminal'

    return featurecloud_app