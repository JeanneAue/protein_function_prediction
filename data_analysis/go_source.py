import os
import requests

def download_file(url, filename):
    """Downloads a file from a URL and saves it to the specified filename."""
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an exception for bad status codes

        with open(filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"Downloaded {filename} successfully.")
    except requests.exceptions.RequestException as e:
        print(f"Error downloading {url}: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage:
url = "http://purl.obolibrary.org/obo/go/go-basic.obo" # Corrected URL
filename = "go-basic.obo"

if not os.path.exists(filename):
    download_file(url, filename)
else:
    print(f"{filename} already exists.")