import json
import requests
from bs4 import BeautifulSoup
import sys
import os
import os.path
import shutil

try : # Tries to create a folder called images
    os.mkdir('images')

except: # If the image folder already exists it deletes the folder and creates a new one to get rid of residual images
    print("Image Folder Already Exists. New Folder Being Made.\n")
    shutil.rmtree('./images') # Deletes the folder
    os.mkdir('images') # Creates a new folder

private_tester = 0 # Create a variable to test if the profile is private or has no posts

# Downloads the image with the filename being the number image it is
def download_image(url, number): # Creates the function and sets its parameters
    image_browser = requests.get(url) # Creates a requests object to go to the image url
    image_file = open("./images/" + str(number) + '.jpg', "w") # Creates a jpg images in the images folder
    image_file.write(image_browser.content) # Writes the content to the image by scraping the image data from the url
    image_file.close() # Closes the file


# Finds the Image links
while private_tester == 0: # Checks to make sure a valid user is entered

    completed = 0 # Value to check when a valid instagram username is entered

    # If it's an invalid instagram username, it prompts you to enter a valid one
    while completed != 200: # While invalid profiles are entered
        link = 'https://www.instagram.com/' + raw_input('Enter Instagram Username: ') # Gets the Instagram profile link
        pagefile = requests.get(link) # Creates an object for the instagram user's information
        completed = pagefile.status_code # Checks if a valid Instagram link was entered. If it was, value should be 200.

        if completed != 200: # If it's an invalid username
            print("\nError: Invalid Instagram Username. Please Enter a New User.\n") # Print this text then prompt for new username

    parser = BeautifulSoup(pagefile.text, 'lxml') # Parses the instagram website using lxml

    ig_json = parser.find('script', text=lambda t: t.startswith('window._sharedData')).text.split(' = ', 1)[1].rstrip(';') # Strips the useless information and just leaves the image data and stores it in a JSON for easy manipulation
    url_file = json.loads(ig_json) # Opens the JSON file and creates an object that we can find the urls from

    for count in url_file['entry_data']['ProfilePage'][0]['graphql']['user']['edge_owner_to_timeline_media']['edges']: # Goes through each instagram post and finds the link to the 480 x 600 image
        photo_link = count['node']['thumbnail_resources'][3]['src'] # Creates a string with the image link
        private_tester += 1 # If it has a post, increment it
        download_image(photo_link, private_tester) # Passes the link and number to the image download function

    if private_tester == 0:
        print("Error: User Private Or Has No Photos. Please Enter A New User.\n") # Prints out error message

print("Images Successfully Downloaded.\n") # Tells us that we downloaded a valid user's urls
