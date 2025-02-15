import os
from Bio import Entrez

# Dynamic DATA_DIR: Creates a "data" folder in your project root directory
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
os.makedirs(DATA_DIR, exist_ok=True)

## you may modify the code to add your own NCBI API KEY !!!   
## you may create a .env file as well 


Entrez.email = "tahagill99@gmail.com"  # this is for users ease , you may put yours

# If you want to use your own email address please uncomment the code below, and put your email in the input of your terminal.
# Entrez.email = input("Enter your NCBI email: ")

            ############################### PLEASE BE CAREFUL OVER HERE [SEE THE COMMENTS] ######################################