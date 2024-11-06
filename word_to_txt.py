# Import the required library
import pywhatkit as pw

# Define the text to be converted to handwriting
txt = """hi, this side mangal chand"""

# Try to generate the handwriting image
try:
    # Convert text to handwriting and save as an image
    pw.text_to_handwriting(txt, "hand_writing.png", [0, 0, 138])
    print("Handwriting image generated successfully.")
except Exception as e:
    # Handle any exceptions that occur
    print(f"An error occurred: {str(e)}")
    print("Consider using alternative libraries or online services.")

# Indicate the end of the script
print("END")