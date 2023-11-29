# MyAppFour - A Template For Federated Learning Application 

Welcome to the appfour template repository! This is designed to help you get started with developing your own federated learning (FL) applications. The template provides a foundational structure, so you can focus on the unique aspects of your application. This README will guide you through using the template, testing your application, and publishing it to the app store.

For registering and testing your apps or using other apps, please visit
[FeatureCloud.ai](https://featurecloud.ai/). And for more information about FeatureCloud architecture,
please refer to 
[The FeatureCloud AI Store for Federated Learning in Biomedicine and Beyond](https://arxiv.org/abs/2105.05734) [[1]](#1).

## Getting Started with Development
1. Create your app
```bash
featurecloud app new YOUR_APPLICATION --template-name app-four
```
Or start by cloning this repository to your local machine.
```bash
git clonehttps://github.com/FeatureCloud/app-four.git
```
Then one can open the p[roject with its favorite IDE and start developing.
2. Understand the Structure:

   MyAppFour class will go through four states for federated scenario and two states for centralized scenario. All FeatureCloud apps 
   have a terminal state which is designed to flag the end of app execution with no other functionality, so, we will not
   consider it as state here. For each state there is counterpart method in MyAppFour which carries all the operations.
   MyAppFour inherits AppFour class, so it would be helpful developers familiarize themselves with `myapp.py` (the primary application logic)
   and `utils.py` (Utility functions and classes to support the app.)

3. Develop Your App:

   We encourage you to follow the following steps to speed up the app development and testing process:
   - Start with centralized scenario in native mode. Test and validate your app in centralized mode.
   - Build the image and test it in containerized mode.
   - continue development by implementing federated simulation scenario in centralized mode.
   - Test you app in centralized mode. From this point, most probably your app will work in the upcoming tests.
   - Optional: Test your app in simulation scenario in containerized mode.
   - Test your app in a real world federated scenario in containerized mode

   Replace the example code in `myapp.py` with your application logic. Ensure you do not rename key classes or methods as other parts of the template depend on them.
   When developing or modifying the application, ensure that the methods of `MyApp` class align with the expected behavior defined in its parent class from `utils.py`.

      Always check the app's configuration and structure against the `config.yml` for consistency.

4. Additional Notes:
   - The application name defined by APP_NAME is case-sensitive. It should contain the name of the app as a string and match the name of the config file.
   - Always name the configuration file as `config.yml`.
   - For Native Centralized mode, data_dir should point to a directory containing all files necessary for centralized training. 
   - For Native Simulation mode, [Details are missing in the provided code]
   - Do not rename the MyApp class and its five methods in myapp.py.
   - For the native execution, if centralized, it points data_dir to a directory containing all files for centralized training.   - 
   - In federated scenarios, for each client an instance of the `MyAppFour` will be executed.
   - The firsty client will be coordinator and the rest will be participants.
   - All the clients will follow the same states from initial to terminal except for the coordinator (first client) who goes to the aggregation state.
   - Data communication happens according to the return of the methods.
   - In `intial` method, What coordinator returns  will be broadcast to all clients. What other clients return will be discarded!
   - Clients receive the broadcast data as `local_training` method's input.
   - As long as `last_round` attribute is `False`, each keep executing `local_training` method for each round.
   - As long as `last_round` attribute is `False`, the coordinator executes `local_training` and `global_aggregation`, respectively will continue .


## How it works
In general, app four supports three scenarios (Centralized, Simulation, and Federated) in two modes (Native and Containerized)
by automatically instantiating and running MyApp. Developers should implement their applications in MyApp class 
which has five methods; except for `centralized` which is dedicated for centralized training, all other predefined
methods will be called in corresponding predefined states.

### Federated
   For federated learning, the app executes the following states:
   1. initial: App runs the `MyApp.initial` method, and broadcasts the returned values of the coordinator to all clients.
   2. local training: App runs the `MyApp.local_training` method, and sends the returned values to coordinator.
   3. global aggregation: App runs the `MyApp.global_aggregation` method, and broadcasts the returned values to clients.
   4. write results: App runs the `MyApp.write_results` method.
   5. terminal: App execution finishes

   ![Federated states diagram](images/federated_state_diagram.png)

   As the federated diagram shows, the app will  start with running initial method/state, and finishes by running write results method/state. 
   Those are the states that will be ran only once. On the other hand, there is a cycle of repetitive transitions between local training and global aggregation
   which goes on until the last round. Each time that the cycle take place will be called a communication round. The cycle can be stopped with
   a fix maximum nuber of rounds and/or stoppage criteria. It is upon developers to stop the cycle via setting `last_round` attribute to `True`.
   
   

### Centralized training
For centralized training, the app will imediately transition from `initial` state to `Centralized`, run the MyApp.centralized
method, and then transition to terminal.

![Centralized states diagram](images/centralized_state_diagram.png)
Beware that `initial` method will not run in centralized scenario. In case you need, you should `initial` method in `centralized`. 
### Modes
FeatureCloud application are primarily designed and implemented to be used in federated workflows on real-world scenarios.
However, for app developers, it would be convenient to implemenyt and test their applications natively, and without containerization.
Therefore, App-four supports tow execution modes to facilitate the app development by supporting most of the functionalities natively.
We suggest, app developers, first develop and test their applications natively, and then test it in containerized mode. 
The modes are completely transparent and one can run the app with terminal without building the docker images. On the other hand, 
they can run their apps using testbed and workflow in containerized mode. 

#### Native
All codes will be run natively, simply by running `python3 main.py` In this mode, no SMPC is supported and required communications
will be handled internally without passing any data on network.
#### Containerized
The app will be executed through the FeatureCloud controller and all the communications passes the controller.
Every app instance will be executed as docker container, where has its own mounted drives. All the functionalities are supported.   

### Scenarios
App-four template supports three scenarios to cover all needs of app developers:
* Centralized: to implement and test centralized model in both modes. Many of the implemented methods, classes, etc., can be used
for both centralized and federated.
* Simulation : it is a federated scenario that can be run in both modes. For containerized, there will be only 
one container instance (one client) that simulates and runs all the clients.
* Federated: Real-world scenario with completely independent app container instances for clients.

All the scenarios are covered in the config file; once there are no centralized or simulation in the config file,
federated scenario will be executed.

## Testing Your Application
Even though FeatureCloud apps are primarily designed and intended to be used in a real-world federated workflow in containerized mode by data owners like hospitals,
after developing you applications you have a couple of options to test it. Here we cover the options based on a recommended order for development (mode vs scenario):

1. Native Centralized:
```bash
$ python PATH/TO/YOR/APP/main.py
```
2. Containerized Centralized:
3. Native Simulation:
```bash
$ python PATH/TO/YOR/APP/main.py
```
4. Containerized Simulation:
5. Containerized federated

All the aforementioned options are for testing you app as a standalone app. Once the app is development is over, to make sure the app 
can be used in workflow (federated project), developers should test their apps alongside other companion apps in a workflow.

## App development
One should start the development with implementing MyApp class, which is inherits FeatureCloudApp.

```python
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

```

   In the default `MyAppFour` class, to break the cycle of rounds (transitions between local training and global aggregation states)
   a `last_round` boolean value was used in a way that coordinator informs others about the last round in  `initial` method,
   while all participants check for last round as their received global parameters in `local_training` method. The value is False, so
   `self.last_round` will be set as `False`. Beware that app-four template always check `self.last_round` after executing the `local_training` 
   method. So, the coordinator transitions to `global_aggregation` state/method, and ignores received `local_parameters` in this case.
   to ensure this is the last round, the coordinator sets `self.last_round` as `True` and broadcasts it. Meantime, participants loop back to 
   `local_training` and await for receiving the global parameters. The moment parameters are received, they run the method and set the last round.
   This time, last round is `True` for all clients (Coordinator and participants), therefore, the cycle is stopped, and all clients will move
   to `write_results` state/method.
   
### Methods
#### Federated
Four methods of `load_data`, `local_training(self, global_parameters)`, `self.last_round = global_parameters`,
`global_aggregation(self, local_parameters)`, and `def write_results` should be implemented.
##### Training loop
the app starts the training loop by transitioning to locat_training state. Then, for every communication round,
coordinator (one of the clients), will transition to aggregate state, runs the aggregation method and broadcasts 
the global parameters to clients. Meantime, other clients, will wait for global parameters in the local training state.
All the transitions and data communication happens automatically. However, developers should expect to receive the global 
parameters and local models as input for corresponding methods.

To break the loop and move to `write_results` method, developers should turn the `last_round` attribute to `True`. 

### Data access
Data will be in different locations, therefore, users should give the path to folder containing data in the config file for centralized and simulation scenario.
For federated scenario, data will be uploaded in a workflow, or path will be provided in the testbed. However, the data root path is different for Native and Containerized mode. 
Developers should use `self.input_root_dir` and `self.output_root_dir` to access the data. 

```python
path_to_train = f"{self.input_root_dir}/{self.config['local_dataset']['train']}"
```
This example not only shows how to access the data, but also shows how one can access the content of config.yml file inside the app.

In a similar fashion, one can handle the output data with this consideration that the app execution in native mode
overrides the output directory.
#### Native mode
##### Simulation data structure
   After running your containers in FeatureCloud test-bed or workflow, the input data will be cloned from the specified path to 
   `/mnt/input` inside the container. However, in native mode the data should be in your app's root directory. Developers can set the name 
   of the directory and subdirectories in the config file:
```yaml
  simulation:
    clients: Client_1, Client_2
    dir: sample_data
    clients_dir: c1,c2 # Comma separated list of clients' directories
```
   Beware that the structure and keys of the simulation configs should not be changed. However, the directory names can be changes. usually,
   in a federated workflow, each clients has its own data, possibly, with different names. Here, for simulation, we assume the 
   data files have the same name but they are provided in different subdirectories. All data are located in `sample_data`, while, for two clients, each client data is located in `c1`
   and `c2`. Beware that for simulation, the number of clients will be according to the number of comma-separated clients ids
   in `clients`, e.g., `Client_1, Client_2`.

##### Centralized input data structure
   The input data structure for centralized scenario in native mode is very similar to the simulation  scenario. Just make
   sure you only put either `centralized` or `simulation` settings in the config file.
```yaml
  centralized:
    data_dir: sample_data/c1
```
Here, we simply used the path to one of the clients' data as input directory.
##### Output data structure
   After running your containers in FeatureCloud test-bed, the results will be cloned from `/mnt/output` in container into the 
   `PATH/TO/CONTROLLER/data/tests` but in native mode, the results wil be written into `PATH/TO/YOUR/APP/results` directory.
   beware that the path is fixed for native mode and the directory will be overwritten after each execution.

#### Containerized mode
   For running a federated workflow or test bed, one container will be instantiated for each client. and the data will be cloned to `mnt/input`.
   However, simulation and centralized scenario in containerized mode, only one container will be ran. So, one should do the folloowinbgs:
   * Run the test-bed for one client
   * Give the path to directory including all clients' data
   * Put config file in the same directory as data.
   * The rest will be handles through the config file similar to native mode. 
     * Put all data in the controller's data directory.
     * Put config.yml file in generic directory.
     * Provide the path to generic directory.

### building the app docker image
Once app implementation is done, building the docker image for testing or adding it to
[FeatureCloud AI store](https://featurecloud.ai/ai-store?view=store&q=&r=0),
developers should provide the following files.
#### Dockerization files

For dockerizing apps, regardless of their applications, there should be some specific files:

1. [Dockerfile](Dockerfile)
2. [server-config](server_config)
   - [docker-entrypoint.sh](server_config/docker-entrypoint.sh)
   - [nginx](server_config/nginx)
   - [supervisord.conf](server_config/supervisord.conf)

Developers should ensure that these files with the same structure and content exist in the same directory as their app
implementation. 


#### App-specific files
All app-specific files should include data or codes strictly dependent on the app's functionality.

##### main.py
Each app should be implemented in a directory that includes the [`main.py`](main.py) file, which in turn comprises either direct
implementation of states or importing them. Moreover, `main` should import `bottle` and `api` packages:
```angular2html
from bottle import Bottle

from api.http_ctrl import api_server
from api.http_web import web_server

import apps.examples.dice

from engine.app import app

server = Bottle()
```
One can implement desired states in [`states.py`](myapp.py) and import it, which because of putting 
[`app_state`](https://github.com/FeatureCloud/FeatureCloud/tree/master/FeatureCloud/app/engine/README.md#registering-states-to-the-app-app_state) on top of state classes, 
merely importing the states and registering them into the [`app` instance](https://github.com/FeatureCloud/FeatureCloud/tree/master/FeatureCloud/app/engine/README.md#app-instance).     

For running the app, inside a docker container, [`app.register()`](https://github.com/FeatureCloud/FeatureCloud/tree/master/FeatureCloud/app/engine/README.md#registering-all-transitions-appregister)
should be called to register and verify all transitions; next, api and servers should mount at corresponding paths; and finally
the server is ready to run the app.

```angular2html
    app.register()
    server.mount('/api', api_server)
    server.mount('/web', web_server)
    server.run(host='localhost', port=5000)
```

All of the codes above, except for importing the app or, alternatively, implementing states, can be exactly same for all apps.  

##### requirements.txt
for installing required python libraries inside the docker image, developers should provide a list of libraries in [requirements.txt](requirements.txt).
Some requirements are necessary for the FeatureCloud library, which should always be listed, are:
```angular2html
bottle
jsonpickle
joblib
numpy
bios
pydot
pyyaml
```

And the rest should be all other app-required libraries.

##### config.yml
Each app may need some hyper-parameters or arguments that the end-users should provide. Such data should be included
in [`config.yml`](https://github.com/FeatureCloud/FeatureCloud/tree/master/FeatureCloud/app#config-file-configyml), which should be read and interpreted by the app. 

### Run YOUR_APPLICATION
## Limitations
   AppFour is limited to have four fixed states for federated scenarion and two fixed states for centralized scenarios.
   This fixed number of states with fixed transitions and data communications causes some limitations:
   * FeatureCloud app plot-states command does not work for apps developed using app-four template.
   * 

#### Prerequisite

To run YOUR_APPLICATION, you should install Docker and FeatureCloud pip package:

```shell
pip install featurecloud
```

Then either download YOUR_APPLICATION image from the FeatureCloud docker repository:

```shell
featurecloud app download featurecloud.ai/YOUR_APPLICATION
```

Or build the app locally:

```shell
featurecloud app build featurecloud.ai/YOUR_APPLICATION
```

Please provide example data so others can run YOUR_APPLICATION with the desired settings in the `config.yml` file.

#### Run YOUR_APPLICATION in the test-bed

You can run YOUR_APPLICATION as a standalone app in the [FeatureCloud test-bed](https://featurecloud.ai/development/test) or [FeatureCloud Workflow](https://featurecloud.ai/projects). You can also run the app using CLI:

```shell
featurecloud test start --app-image featurecloud.ai/YOUR_APPLICATION --client-dirs './sample/c1,./sample/c2' --generic-dir './sample/generic'
```

### Future Plans:

* Extend the state transition mechanism for more complex training scenarios.
* Integrate with additional machine learning frameworks.

### References
<a id="1">[1]</a> 
Matschinske, J., Späth, J., Nasirigerdeh, R., Torkzadehmahani, R., Hartebrodt, A., Orbán, B., Fejér, S., Zolotareva,
O., Bakhtiari, M., Bihari, B. and Bloice, M., 2021.
The FeatureCloud AI Store for Federated Learning in Biomedicine and Beyond. arXiv preprint arXiv:2105.05734.
