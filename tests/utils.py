
import sys




def print_col(text, color):

    if(color == "black"):

        print("\033[30m{}\033[0m".format(text))

    elif(color == "red"):

        print("\033[31m{}\033[0m".format(text))

    elif(color == "green"):

        print("\033[32m{}\033[0m".format(text))

    elif(color == "yellow"):

        print("\033[33m{}\033[0m".format(text))

    elif(color == "blue"):

        print("\033[34m{}\033[0m".format(text))

    elif(color == "magenta"):

        print("\033[35m{}\033[0m".format(text))

    elif(color == "cyan"):

        print("\033[36m{}\033[0m".format(text))

    elif(color == "white"):

        print("\033[37m{}\033[0m".format(text))

    else:

        print("{}".format(text))


# ---

def stdout_col(text, color):

    if(color == "black"):

        sys.stdout.write("\033[30m{}\033[0m".format(text))

    elif(color == "red"):

        sys.stdout.write("\033[31m{}\033[0m".format(text))

    elif(color == "green"):

        sys.stdout.write("\033[32m{}\033[0m".format(text))

    elif(color == "yellow"):

        sys.stdout.write("\033[33m{}\033[0m".format(text))

    elif(color == "blue"):

        sys.stdout.write("\033[34m{}\033[0m".format(text))

    elif(color == "magenta"):

        sys.stdout.write("\033[35m{}\033[0m".format(text))

    elif(color == "cyan"):

        sys.stdout.write("\033[36m{}\033[0m".format(text))

    elif(color == "white"):

        sys.stdout.write("\033[37m{}\033[0m".format(text))

    else:

        sys.stdout.write("{}".format(text))

# ---


