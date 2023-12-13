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
import os
import shutil
import bios
import ast
import abc
from time import sleep
from FeatureCloud.app.engine.app import App, app
import states
from bottle import Bottle
from FeatureCloud.app.api.http_ctrl import api_server
from FeatureCloud.app.api.http_web import web_server
import yaml


def print_configurations(data):
    print("\t\t\tConfigurations...")
    print("\n****************************************\n")
    formatted_yaml = yaml.dump(data, indent=2)
    print(formatted_yaml)
    print("\n****************************************\n")


def is_native():
    '''The function checks if the "PATH_PREFIX" environment variable is set and returns True if it is not
    set, indicating that the code is running in a native environment.
    
    Returns
    -------
        a boolean value. If the `path_prefix` variable is not empty, it will return `False`. Otherwise, it
    will return `True`.
    
    '''
    path_prefix = os.getenv("PATH_PREFIX")
    if path_prefix:
        return False
    return True


CONFIG_FILE = "config.yml"


def read_config():
    '''The function `read_config()` reads a configuration file either from the native file system or from a
    specific directory if it is not native.
    
    Returns
    -------
        the contents of the config file.
    
    '''
    file_path = CONFIG_FILE
    if not is_native():
        file_path = f"app/{CONFIG_FILE}"
    return bios.read(file_path)


def is_simulation(app_name):
    '''The function `is_simulation` checks if a given app name is configured for simulation.
    
    Parameters
    ----------
    app_name
        The `app_name` parameter is a string that represents the name of the application.
    
    Returns
    -------
        a boolean value. It returns True if the 'simulation' key exists in the configuration for the given
    app_name and has a non-null value. Otherwise, it returns False.
    
    '''
    config = read_config()[app_name]
    return config.get('simulation', None) is not None


def get_root_dir(input_dir=True, simulation_dir=None):
    '''The function `get_root_dir` returns the root directory path based on the input directory flag and
    simulation directory.
    
    Parameters
    ----------
    input_dir, optional
        A boolean flag indicating whether the function should return the input directory or the output
    directory. If `input_dir` is `True`, the function will return the input directory. If `input_dir` is
    `False`, the function will return the output directory.
    simulation_dir
        The `simulation_dir` parameter is a string that represents the directory path within the root
    directory. It is an optional parameter and can be set to `None` if not needed.
    
    Returns
    -------
        a string that represents the root directory path.
    
    '''
    simulation_dir = "/" + simulation_dir if simulation_dir else ''

    if input_dir:
        if is_native():
            return f".{simulation_dir}"
        return f"/mnt/input{simulation_dir}"
    if is_native():
        return f"./results{simulation_dir}"
    return f"/mnt/output{simulation_dir}"


def is_centralized(app_name: str):
    '''The function `is_centralized` checks if a given application name is configured as centralized in a
    configuration file.
    
    Parameters
    ----------
    app_name : str
        A string representing the name of the application.
    
    Returns
    -------
        the value of the 'centralized' key in the config dictionary for the given app_name. If the key is
    not present, it will return False.
    
    '''
    config = read_config()[app_name]
    return config.get('centralized', False)


def file_has_class(file_path: str, class_name: str):
    '''The function `file_has_class` checks if a given file contains a class with a specific name.
    
    Parameters
    ----------
    file_path : str
        The `file_path` parameter is a string that represents the path to the file that you want to check
    for the presence of a specific class.
    class_name : str
        The `class_name` parameter is a string that represents the name of the class you want to check for
    in the file.
    
    Returns
    -------
        a boolean value. It returns True if the specified file contains a class with the given class name,
    and False otherwise.
    
    '''
    with open(file_path, 'r') as file:
        source_code = file.read()

    try:
        parsed = ast.parse(source_code)
    except SyntaxError:
        return False

    for node in ast.walk(parsed):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            return True

    return False


# The class `AppFour` is an abstract base class that defines the structure and methods for a
# simulation application with a centralized coordinator and multiple clients.
class AppFour(abc.ABC):
    def __init__(self, config, simulation_dir=None):
        self.config = config
        self.input_root_dir = get_root_dir(simulation_dir=simulation_dir)
        self.output_root_dir = get_root_dir(input_dir=False, simulation_dir=simulation_dir)
        if os.path.exists(self.output_root_dir):
            shutil.rmtree(self.output_root_dir)
        os.makedirs(self.output_root_dir)

        if 'logic' in self.config:
            self.mode = self.config['logic']['mode']
            self.max_na_rate = self.config['logic']['max_na_rate']

        self.log(f"Input directory: {self.input_root_dir},\n",
                 f"Output directory: {self.output_root_dir}\n",
                 f"Mode: {self.mode}, Max NA rate: {self.max_na_rate}")

        self.last_round = False

    @abc.abstractmethod
    def initial(self):
        """

        Returns
        -------
        data_to_broadcast: list
#             Values or messages from coordinator to all client
        """
        pass

    @abc.abstractmethod
    def local_training(self, global_parameters: list):
        """

        Parameters
        ----------
        global_parameters

        Returns
        -------

        """
        pass

    @abc.abstractmethod
    def global_aggregation(self, local_parameters: list):
        """

        Parameters
        ----------
        local_parameters

        Returns
        -------

        """
        pass

    @abc.abstractmethod
    def write_results(self):
        """

        Returns
        -------

        """
        pass

    @abc.abstractmethod
    def centralized(self):
        pass


def run_client(app: app, client_id: str, clients_ids:list, coordinator_role: bool, shared_memory: dict):
    '''The function `run_client` registers the client with the app, handles the setup, and assigns the
    shared memory.
    
    Parameters
    ----------
    app : app
        The "app" parameter is an instance of a class that represents the application you are running. It
    likely has methods and attributes that allow you to interact with the application.
    client_id : str
        The client_id parameter is a string that represents the unique identifier for the client.
    clients_ids : list
        The `clients_ids` parameter is a list of client IDs. It represents the IDs of all the clients that
    are participating in the application.
    coordinator_role : bool
        The `coordinator_role` parameter is a boolean value that indicates whether the client is acting as
    the coordinator or not. If `coordinator_role` is `True`, it means the client is the coordinator. If
    it is `False`, it means the client is not the coordinator.
    shared_memory : dict
        The `shared_memory` parameter is a dictionary that is used to share data between different clients
    in the application. It allows clients to read and write data that can be accessed by other clients.
    
    '''
    app.register()
    app.handle_setup(client_id=client_id,
                     coordinator=coordinator_role,
                     clients=clients_ids)
    app.shared_memory = shared_memory


# The `Controller` class is responsible for managing and coordinating communication between multiple
# clients and an application.
class Controller:
    def __init__(self, clients_id: str):
        '''The function initializes a dictionary with client IDs as keys and empty lists as values.
        
        Parameters
        ----------
        clients_id : str
            A string representing the IDs of the clients.
        
        '''
        self.clients = {id.strip(): [] for id in clients_id}

    def register(self, client: str, app: app, coordinator: bool):
        '''The function registers a client with an app and updates the list of clients.
        
        Parameters
        ----------
        client : str
            The `client` parameter is a string that represents the client's name or identifier.
        app : app
            The "app" parameter is an instance of a class or object that represents an application. It is
        used to perform various operations related to the application, such as registering and setting
        up the client.
        coordinator : bool
            The "coordinator" parameter is a boolean value that indicates whether the client being
        registered is a coordinator or not.
        
        '''
        app.register()
        app.handle_setup(client_id=client,
                         coordinator=coordinator,
                         clients=list(self.clients.keys()))
        self.clients[client] = app

    def run(self):
        '''The function runs a loop that checks for available data from clients and sends it to the
        appropriate destination client or coordinator.
        
        '''
        while True:
            for client in self.clients:
                if self.data_available(client):
                    data, dest_client = self.check_outbound(client)
                    if dest_client:
                        print("send_to_participant")
                        self.set_inbound(data, source_client=client, dest_client=dest_client)
                    else:
                        if self.clients[client].coordinator:
                            # broadcast
                            print("broadcast", list(self.clients.keys())[1:])
                            for dest_client in list(self.clients.keys())[1:]:
                                self.set_inbound(data, source_client=client, dest_client=dest_client)
                        else:
                            print("send_to_coordinator")
                            self.set_inbound(data, source_client=client, dest_client=list(self.clients.keys())[0])

            sleep(1)
            if self.finished():
                break

    def check_outbound(self, client: str):
        '''The function "check_outbound" returns the outgoing data and destination for a given client.
        
        Parameters
        ----------
        client : str
            The client parameter is a string that represents the name or identifier of the client for which
        we want to check the outbound status.
        
        Returns
        -------
            two values: `data` and `dest`.
        
        '''
        dest = self.status(client)['destination']

        data = self.clients[client].handle_outgoing()
        return data, dest

    def set_inbound(self, data: dict or list, source_client: str, dest_client: str):
        '''The function sets the inbound data from a source client to a destination client.
        
        Parameters
        ----------
        data : dict or list
            The `data` parameter is either a dictionary or a list. It represents the incoming data that
        needs to be processed by the `handle_incoming` method of the `dest_client`.
        source_client : str
            The `source_client` parameter is a string that represents the client from which the data is
        being sent.
        dest_client : str
            The `dest_client` parameter is a string that represents the destination client. It is used to
        identify the client that will handle the incoming data.
        
        '''
        self.clients[dest_client].handle_incoming(data, source_client)

    def finished(self):
        '''The function checks if all the clients have finished their tasks.
        
        Returns
        -------
            a boolean value indicating whether all the clients' status is marked as "finished".
        
        '''
        finished = [self.status(app)['finished'] for app in self.clients]
        return all(finished)

    def status(self, client: str):
        '''The function `status` returns the status of a client in a dictionary format.
        
        Parameters
        ----------
        client : str
            The `client` parameter is a string that represents the name or identifier of a client. It is
        used to retrieve the status of the specified client from the `self.clients` dictionary.
        
        Returns
        -------
            a dictionary containing various status information about the specified client.
        
        '''
        app = self.clients[client]
        client_status = {
            'available': app.status_available,
            'finished': app.status_finished,
            'message': app.status_message if app.status_message else (
                app.current_state.name if app.current_state else None),
            'progress': app.status_progress,
            'state': app.status_state,
            'destination': app.status_destination,
            'smpc': app.status_smpc,
        }
        return client_status

    def data_available(self, client: str):
        '''The function "data_available" returns the availability status of a client.
        
        Parameters
        ----------
        client : str
            The "client" parameter is a string that represents the name or identifier of a client.
        
        Returns
        -------
            the value of the 'available' key from the status of the client.
        
        '''
        return self.status(client)['available']


def simulate(config: dict, app_name: str, MyApp: AppFour):
    '''The `simulate` function takes in a configuration dictionary, an app name, and a `MyApp` object, and
    simulates the execution of a feature cloud application using the provided configuration.
    
    Parameters
    ----------
    config : dict
        A dictionary containing the configuration settings for the simulation.
    app_name : str
        The `app_name` parameter is a string that represents the name of the application being simulated.
    MyApp : AppFour
        `MyApp` is a class that represents the feature cloud application. It is used to create an instance
    of the application with the given configuration and simulation directory.
    
    '''
    print(f"Simulated: {config['simulation']}")
    clients_ids = [id.strip() for id in config['simulation']['clients'].split(',')]
    clients_dirs = config['simulation']['clients_dir'].split(',')
    simulation_dir = config['simulation'].get('dir', None)
    if simulation_dir:
        clients_dirs = [f"{simulation_dir}/{d}" for d in clients_dirs]
    controller = Controller(clients_ids)
    for i, (client_id, client_dir) in enumerate(zip(clients_ids, clients_dirs)):
        app = App()
        app_four = MyApp(config=read_config()[app_name], simulation_dir=client_dir)
        states.create_client(app_four, app)
        controller.register(client_id, app, coordinator=i == 0)
    controller.run()


def centralized(app_class: AppFour, app_name: str):
    '''The function `centralized` creates a centralized instance of a feature cloud app and registers it
    with a general app instance.
    
    Parameters
    ----------
    app_class : AppFour
        The app_class parameter is the class of the feature cloud app that you want to create an instance
    of. It should be a subclass of FeatureCloudApp.
    app_name : str
        The `app_name` parameter is a string that represents the name of the application.
    
    '''
    config = read_config()[app_name]
    print(f"centralized: {config['centralized']}")
    featurecloud_app = App()
    app_four = app_class(config=config, simulation_dir=config['centralized'].get('data_dir', None))
    featurecloud_app = states.create_client(app_four, featurecloud_app, centralized=True)
    featurecloud_app.register()
    featurecloud_app.handle_setup(client_id='1', coordinator=True, clients=['1'])


def federated(app_class: AppFour, app_name: str):
    '''The `federated` function sets up a federated learning server for a given app.
    
    Parameters
    ----------
    app_class : AppFour
        The app_class parameter is the class that represents the feature cloud app. It should be a subclass
    of FeatureCloudApp.
    app_name : str
        The `app_name` parameter is a string that represents the name of the application. It is used to
    retrieve the configuration for the specified application from the `read_config()` function.
    
    '''
    my_app = app_class(read_config()[app_name])
    states.generate_federated_states(app_four=my_app)
    app.register()
    server = Bottle()
    server.mount('/api', api_server)
    server.mount('/web', web_server)
    server.run(host='localhost', port=5000)