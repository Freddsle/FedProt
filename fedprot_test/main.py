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
import utils
from myapp import MyAppFour, APP_NAME

if __name__ == '__main__':
    config = utils.read_config()[APP_NAME]
    utils.print_configurations(config)
    if utils.is_native():
        print("The app will run in Native mode...")
        if utils.is_centralized(APP_NAME):
            print("Centralized analysis...")
            utils.centralized(MyAppFour, APP_NAME)
        elif utils.is_simulation(APP_NAME):
            print("Simulating federated analysis...")
            utils.simulate(config, APP_NAME, MyAppFour)
        else:
            raise NotImplemented(f"Native execution is only available for 'centralized' or `simulation` scenarios")
    else:
        if utils.is_centralized(APP_NAME):
            utils.centralized(MyAppFour, APP_NAME)
        elif utils.is_simulation(APP_NAME):
            utils.simulate(config, APP_NAME, MyAppFour)
        else:
            utils.federated(MyAppFour, APP_NAME)
