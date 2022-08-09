8  # Helper functions for app.py

from itertools import chain
from flask import render_template, request


def error(message, code=400):
    """Render error message"""
    return render_template("error.html", message=message), code


def digest(seq, rec1, rec2):
    frags = digest_first_enzyme(rec1, seq)
    if rec2:
        digest_second_enzyme(frags, rec2)
        frags = list(chain.from_iterable(frags))

    # Connect beginning and end of sequence if circular topology
    topology = request.form.get("topology-options")
    if topology == "circular":
        frags[len(frags) - 1] = frags[len(frags) - 1] + frags[0]
        del frags[0]

    sort_order = sort_fragments(frags)
    frag_info = get_frag_info(frags)

    return frag_info, sort_order


def digest_first_enzyme(rec1, seq):
    # Cut target sequence into fragments with first recognition sequence
    frags = seq.split(rec1.replace("*", ""))
    # Since split() removed recognition sequence, add it back to each fragment
    rec1_5_prime = rec1.split("*")[0]
    rec1_3_prime = rec1.split("*")[1]
    for frag in range(len(frags)):
        if frag == 0:
            frags[frag] = frags[frag] + rec1_5_prime
        elif frag == len(frags) - 1:
            frags[frag] = rec1_3_prime + frags[frag]
        else:
            frags[frag] = rec1_3_prime + frags[frag] + rec1_5_prime
    return frags


def digest_second_enzyme(frags, rec2):
    rec2_5_prime = rec2.split("*")[0]
    rec2_3_prime = rec2.split("*")[1]
    for i in range(len(frags)):
        frags[i] = frags[i].split(rec2.replace("*", ""))

        # Add recognition sequence back to fragments
        if type(frags[i] == "list"):  # Split produces sublist
            for subfrag in range(len(frags[i])):
                if frags[i][subfrag] == 0:
                    frags[i][subfrag] = frags[i][subfrag] + rec2_5_prime
                elif frags[i][subfrag] == len(frags[i]) - 1:
                    frags[i][subfrag] = rec2_5_prime + frags[i][subfrag]
                else:
                    frags[i][subfrag] = (rec2_3_prime
                                         + frags[i][subfrag] + rec2_5_prime)


def sort_fragments(frags):
    sort_order = request.form.get("sort-options")
    if sort_order == "largest-to-smallest":
        frags.sort(reverse=True, key=lambda x: len(x))
    if sort_order == "smallest-to-largest":
        frags.sort(key=lambda x: len(x))
    return sort_order


def get_frag_info(frags):
    frag_info = []
    for i in range(len(frags)):
        frag_info.append(
            {'frag_num': i + 1,
             'seq': frags[i],
             'size': len(frags[i])})
    return frag_info
