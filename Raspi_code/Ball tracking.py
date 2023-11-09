import cv2
import time
import numpy as np
import csv

print(cv2.__version__)

# Set the camera resolution (adjust these values as needed)
width = 640  # Width in pixels
height = 480  # Height in pixels
# Initialize the camera with the specified resolution
cap = cv2.VideoCapture(0)

# set resolution 
cap.set(3, width)  # Set the width
cap.set(4, height)  # Set the height

# set desired frame rate
cap.set(cv2.CAP_PROP_FPS, 60)

lower_green = np.array([40, 100, 40])
upper_green = np.array([80, 255, 255])

min_area = 2000
max_area = 8000

# Check if the camera opened successfully
if not cap.isOpened():
    print("Error: Could not open camera.")
    exit()

frame_count = 0
start_time = time.time()


while True:
    # Capture a frame
    ret, frame = cap.read()

    # Check if the frame was captured successfully
    if not ret:
        print("Error: Could not read frame.")
        break

    hsv_frame = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV)

    mask = cv2.inRange(hsv_frame, lower_green, upper_green)

    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    centroid_x, centroid_y = 0, 0
    for contour in contours:
        area = cv2.contourArea(contour)
        if min_area < area < max_area:
            x, y, w, h = cv2.boundingRect(contour)
            cv2.rectangle(frame, (x, y), (x + w, y + h), (0, 255, 0), 2)
            centroid_x = int((x + x + w) / 2)
            centroid_y = int((y + y + h) / 2)

    cv2.circle(frame, (centroid_x, centroid_y), 5, (0, 0, 255), -1)

    # Display the captured frame
    cv2.imshow('Green Detection', frame)
    print(f"Center X: {centroid_x}, Center Y: {centroid_y}")
    frame_count += 1
    if frame_count % 30 == 0:
        end_time = time.time()
        elapsed_time = end_time - start_time
        frame_rate = frame_count / elapsed_time
        print(f"Frame Rate: {frame_rate:.2f} fps")


    # Break the loop if 'q' key is pressed
    if cv2.waitKey(1) & 0xFF == ord('q'):
        break

# Release the camera and close all OpenCV windows
cap.release()
cv2.destroyAllWindows()
