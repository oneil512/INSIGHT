from colorama import Fore
import os
import json

def get_objective(path):
    state_file_path = os.path.join(path, "state.json")

    # Check if the file exists
    if not os.path.isfile(state_file_path):
        print(f"Error: path {path} does not contain a state.json file.")
        exit()

    # Read the contents of the state.json file
    with open(state_file_path, "r") as state_file:
        try:
            state_data = json.load(state_file)
        except json.JSONDecodeError:
            print("Error: The state.json file is not valid JSON.")
            exit()

    # Check if the objective key exists in the state data
    if "objective" not in state_data:
        print("Error: The objective key is not present in the state.json file.")
        exit()

    return state_data["objective"]
   

def get_input(prompt, type_=None, min_=None, max_=None, range_=None):
    if min_ is not None and max_ is not None and max_ < min_:
        raise ValueError("min_ must be less than or equal to max_.")
    while True:
        ui = input(prompt)
        if type_ is not None:
            try:
                ui = type_(ui)
            except ValueError:
                print(f"Input type must be {type_.__name__}!")
                continue
        if max_ is not None and ui > max_:
            print(f"Input must be less than or equal to {max_}.")
        elif min_ is not None and ui < min_:
            print(f"Input must be greater than or equal to {min_}.")
        elif range_ is not None and ui not in range_:
            if isinstance(range_, range):
                template = "Input must be between {} and {}."
                print(template.format(range_.start, range_.stop))
            else:
                template = "Input must be {}."
                print(template.format(", ".join(map(str, range_))))
        else:
            return ui


def prompt_user():
    reload_path = None
    execution_type = get_input("Would you like to run a new execution (1) or resume an old execution (2)? ", type_=int, range_=(1, 2))

    if execution_type == 1:
        objective = get_input("Enter your objective: ")
        print(f"Your objective is: {objective}")
    elif execution_type == 2:
        while True:
            reload_path = get_input("Enter the path to your previous execution (e.g. out/Cure breast cancer_2023-04-29_15-13-10): ")

            if os.path.isdir(reload_path):
                break
            else:
                print("Directory does not exist. Please try again.")

        objective = get_objective(reload_path)
        print(f"Resuming execution from: {reload_path}")

    print("Now we will do tool selection.")
    tools = ["MYGENE", "PUBMED"]
    tool_flags = {}

    for tool in tools:
        while True:
            tool_prompt = f"Do you want to use {tool}? Type 1 for yes and 0 for no: "
            tool_input = get_input(tool_prompt, type_=int, range_=(0, 1))
            if tool_input is not None:
                tool_flags[tool] = bool(tool_input)
                break
            else:
                print(f"Unrecognized input, defaulting to 'yes' for using {tool}")
                tool_flags[tool] = True

    iterations_prompt = "How many iterations would you like to run? "
    iterations = get_input(iterations_prompt, type_=int, min_=1)

    document_check = get_input("Would you like to load your own data as a document? Type 1 for yes and 0 for no: ", type_=int, range_=(0, 1))
    document_path = None
    if document_check == 1:
        while True:
            document_path = get_input("Enter the path to the document: ")
            if os.path.isfile(document_path):
                break
            else:
                print("File does not exist. Please try again.")

    print("\nHere are the options you've selected:")
    print(f"Objective: {objective}")
    for tool, flag in tool_flags.items():
        print(f"Use {tool}: {'Yes' if flag else 'No'}")
    print(f"Iterations: {iterations}")
    if document_path is not None:
        print(f"Document: {document_path}")

    correct_prompt = "Does this look correct? Type 1 for yes and 0 for no: "
    correct = get_input(correct_prompt, type_=int, range_=(0, 1))

    if correct:
        print(Fore.GREEN + "\033[1mStarting INSIGHT with your options!\033[0m")
    else:
        prompt_user()  # Start over if the user is not satisfied

    return objective, tool_flags, iterations, reload_path, document_path