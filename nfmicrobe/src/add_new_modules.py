#!/usr/bin/env python

import argparse
import fileinput
import json
import sys


def parse_args(args=None):
    description = "Add a new module to nf-microbe"
    epilog = "Example usage: python add_new_module.py --module fastqc"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-d",
        "--directory",
        default=".",
        help="Path to the directory containing nf-microbe/modules repository.",
    )
    parser.add_argument(
        "-m",
        "--module",
        help="Name of module to add",
    )
    return parser.parse_args(args)


# add module to modules.json
def update_modules_json(module_name, modules_json="modules.json"):
    # load json file
    with open(modules_json) as file:
        modules = json.load(file)
    # add module to json
    modules["repos"]["https://github.com/nf-core/modules.git"]["modules"]["nf-core"][module_name] = modules["repos"][
        "https://github.com/nf-core/modules.git"
    ]["modules"]["nf-core"]["coverm/contig"]
    # sort json dict
    sorted_modules = modules.copy()
    sorted_modules["repos"]["https://github.com/nf-core/modules.git"]["modules"]["nf-core"] = dict(
        sorted(modules["repos"]["https://github.com/nf-core/modules.git"]["modules"]["nf-core"].items())
    )
    # write json file
    with open(modules_json, "w") as file:
        json.dump(sorted_modules, file, indent=4)


# add module name to .nf-core.yml to ignore updates
def add_module_to_nfcore_yml(module_name, nfcore_yml=".nf-core.yml"):
    with fileinput.FileInput(nfcore_yml, inplace=True) as file:
        update = False
        repo = False
        nfcore = False
        for line in file:
            if "update:" in line:
                update = True
            if "https://github.com/nf-core/modules.git:" in line:
                repo = True
            if "nf-core:" in line:
                nfcore = True
            if update and repo and nfcore:
                print(line.replace("nf-core:\n", "nf-core:\n      " + module_name + ": False\n"), end="")
                nfcore = False
            else:
                print(line, end="")


# main function
def main(args=None):
    args = parse_args(args)

    update_modules_json(args.module, args.directory + "/modules.json")

    add_module_to_nfcore_yml(args.module, args.directory + "/.nf-core.yml")


if __name__ == "__main__":
    sys.exit(main())
