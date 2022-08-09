# Fledask framework for backend coding of Restriction Digest Simulator website

from flask import Flask, render_template, request
from flask_session import Session
from tempfile import mkdtemp

from enzymes import ENZYMES
from helpers import error, digest

# Configure application
app = Flask(__name__)

# Ensure templates are auto-reloaded
app.config["TEMPLATES_AUTO_RELOAD"] = True


# Ensure responses aren't cached
@app.after_request
def after_request(response):
    response.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
    response.headers["Expires"] = 0
    response.headers["Pragma"] = "no-cache"
    return response


# Configure session to use filesystem (instead of signed cookies)
app.config["SESSION_FILE_DIR"] = mkdtemp()
app.config["SESSION_PERMANENT"] = False
app.config["SESSION_TYPE"] = "filesystem"
Session(app)


@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "GET":
        return render_template("index.html", enzymes=ENZYMES.items())
    else:
        # Get DNA sequence to be digested
        seq_input = request.form.get("seq_input")
        seq = seq_input.upper().replace("\n", "").replace("\r", "")
        for char in seq:
            if char not in ['A', 'T', 'C', 'G']:
                return error(
                    "Entered sequence has one or more invalid characters", 403)

        # Get restriction enzymes or custom recognition sequence
        rec1, rec2 = get_enzymes()
        if not rec1 and not rec2:
            return error("Choose at least one restriction enzyme", 403)

        # Check recognition sequence contained within target sequence
        if rec1.replace("*", "") not in seq or rec2.replace("*", "") not in seq:
            return error("Recognition sequence(s) not found in target sequence")
        if not rec1 and rec2 is not None:
            rec1 = rec2

        frag_info, sort_order = digest(seq, rec1, rec2)

        return render_template(
            "result.html", frag_info=frag_info, sort_order=sort_order)


def get_enzymes():
    rec1 = request.form.get("restriction-enzyme-1")
    if not rec1:
        rec1 = request.form.get("custom-rec-1").upper()
    rec2 = request.form.get("restriction-enzyme-2")
    if not rec2:
        rec2 = request.form.get("custom-rec-2").upper()
    return rec1, rec2
