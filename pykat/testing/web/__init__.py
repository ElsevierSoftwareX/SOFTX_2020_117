from flask import Flask
import os

global app
app = Flask(__name__, instance_path=os.getcwd())