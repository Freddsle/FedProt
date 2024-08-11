import os
import logging
import time
from dataclasses import dataclass
from typing import List

from FeatureCloud.api.imp.test import commands as fc
from FeatureCloud.api.imp.controller import commands as fc_controller

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

@dataclass
class Experiment():
    clients: List[str] # List of client folders (absolute paths)
    app_image_name = "fedprot" # The name of the app image to be used,  
    query_interval = 1  # The query interval in seconds, 1 for default
    controller_host = None # The controller host address, None for default (http://localhost:8000)


def start_controller(data_dir):
    """
    Starts the FeatureCloud controller with the given data directory.
    Args:
        data_dir (str): Absolute path to the directory containing all client files (folders).
    """
    if not os.path.exists(data_dir):
        raise FileNotFoundError(f"Data directory {data_dir} does not exist!")
    
    # Stop the controller if running and start it again
    try:
        fc_controller.stop(name="")
        logging.info("Controller stopped successfully.")
    except Exception as e:
        logging.error(f"Could not stop the controller: {e}")
        raise RuntimeError(f"Could not stop the controller as needed for a clean start: {e}")

    logging.info(f"Starting the FeatureCloud controller with the data directory {data_dir}")
    try:
        fc_controller.start(name="fc-controller", port=8000, data_dir=data_dir, 
                            controller_image="", with_gpu=False, mount="", blockchain_address="")
    except Exception as e:
        logging.error(f"Could not start the controller: {e}")
        raise RuntimeError(f"Could not start the controller! Error: {e}")
    logging.info("Controller started successfully.")
    
    logging.info("Waiting 3s to ensure the controller is properly running...")
    time.sleep(3)

def run_test(experiment, data_dir, retry=5):
    """
    Runs the given experiment.
    Args:
        experiment (Experiment): The experiment to be run.
        data_dir (str): The data directory used by the controller.
        retry (int): Number of retries in case of failure.
    """
    if retry == 0:
        logging.info("_______________EXPERIMENT_______________")
    
    exp_id = start_test(experiment=experiment, data_dir=data_dir, retry=retry)
    instances = check_test(exp_id=exp_id, experiment=experiment, data_dir=data_dir, retry=retry)

    result_files = [f"results_test_{exp_id}_client_{info['id']}_{info['name']}.zip" for info in instances]
    return result_files

def start_test(experiment, data_dir, retry):
    """
    Starts a test with the given experiment configuration.
    Args:
        experiment (Experiment): The experiment to be started.
        data_dir (str): The data directory used by the controller.
        retry (int): The number of retries in case of failure.
    Returns:
        int: The ID of the started test.
    """
    if not experiment.clients:
        raise RuntimeError("The clients should be given!")
    if not experiment.app_image_name:
        raise RuntimeError("The app image name should be given!")
    
    experiment.controller_host = experiment.controller_host or "http://localhost:8000"  
    # experiment.channel = experiment.channel or "local"
    # experiment.query_interval = experiment.query_interval or 5

    client_dirs = ",".join([os.path.relpath(client, data_dir) for client in experiment.clients])
    
    try:
        exp_id = fc.start(
            controller_host=experiment.controller_host,
            app_image=experiment.app_image_name,
            client_dirs=client_dirs,
            channel="local",
            generic_dir="",
            query_interval=1,
            download_results="tests"
        )
        logging.info("Test started successfully!")
    except Exception as e:
        retry += 1
        if retry > 5:
            logging.error(f"Test could not be started after 5 retries: {e}")
            raise RuntimeError(f"Test could not be started more than 5 times! Error: {e}")
        
        logging.warning(f"Could not start the test! Error: {e}. Retrying for the {retry}nth time.")
        time.sleep(5)
        
        start_controller(data_dir=data_dir)
        return start_test(experiment=experiment, data_dir=data_dir, retry=retry)
    
    return int(exp_id)

def check_test(exp_id, experiment, data_dir, retry):
    """
    Checks the status of the test until it is finished.
    Args:
        exp_id (int): The ID of the test to be checked.
        experiment (Experiment): The experiment that was run.
        data_dir (str): The data directory used by the controller.
        retry (int): The number of current retries in case of failure.
    Returns:
        list[dict]: A list of dictionaries containing the instances (clients) of the test.
    """
    experiment.controller_host = experiment.controller_host or "http://localhost:8000"
    
    while True:
        try:
            test_info = fc.info(controller_host=experiment.controller_host, test_id=exp_id)
        except Exception as e:
            logging.error(f"Could not get the test info: {e}")
            retry += 1
            logging.warning(f"Retrying for the {retry}nth time.")
            if retry > 5:
                logging.error(f"Test could not be checked after 5 retries: {e}")
                raise RuntimeError(f"Test could not be checked more than 5 times! Error: {e}")
            time.sleep(5)
            return run_test(experiment=experiment, data_dir=data_dir, retry=retry)
        
        status = test_info.iloc[0]['status']
        instances = test_info.iloc[0]['instances']
        
        if status == "finished":
            logging.info("Test finished successfully!")
            return instances
        elif status in {"error", "stopped"}:
            logging.error(f"Test finished with an error or was stopped! Status: {status}")
            raise RuntimeError(f"Test finished with an error or was stopped! Status: {status}")
