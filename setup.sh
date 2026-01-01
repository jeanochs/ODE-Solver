#!/bin/bash
/usr/bin/python3 -m venv venv_main

# Download the packages
source venv_main/bin/activate
pip install -r requirements_main.txt
