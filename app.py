import dash
from flask import Flask, send_from_directory

import nltk
import ssl

try: 
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    pass
else:
    ssl._create_default_https_context = _create_unverified_https_context

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']


# Normally, Dash creates its own Flask server internally. By creating our own,
# we can create a route for downloading files directly:
server = Flask(__name__)

app = dash.Dash(__name__, external_stylesheets=external_stylesheets,server=server)
app.config.suppress_callback_exceptions = True
print("***** Inside app.py *****")
