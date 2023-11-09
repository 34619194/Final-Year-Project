import cv2
import numpy as np

# Set the camera resolution
width = 640
height = 480

# Initialize the camera with the specified resolution
cap = cv2.VideoCapture(0)
cap.set(3, width)
cap.set(4, height)

lower_green = np.array([40, 100, 40])
upper_green = np.array([80, 255, 255])

# Check if the camera opened successfully
if not cap.isOpened():
    print("Error: Could not open camera.")
    exit()

while True:
    # Capture a frame
    ret, frame = cap.read()

    if not ret:
        print("Error: Could not read frame.")
        break
    original_frame = frame.copy()  # Create a copy of the original frame
    
    hsv_frame = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV)

    mask = cv2.inRange(hsv_frame, lower_green, upper_green)
    
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    for contour in contours:
        x, y, w, h = cv2.boundingRect(contour)
        cv2.rectangle(frame, (x, y), (x + w, y + h), (0, 255, 0), 2)
        
    # Create a side-by-side display
    combined_image = cv2.hconcat([original_frame, frame])

    # Display the combined image
    cv2.imshow('Original vs. Green Mask', combined_image)

    # Break the loop if 'q' key is pressed
    if cv2.waitKey(1) & 0xFF == ord('q'):
        break

# Release the camera and close all OpenCV windows
cap.release()
cv2.destroyAllWindows()
