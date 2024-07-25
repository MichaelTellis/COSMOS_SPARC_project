import qrcode
import matplotlib.pyplot as plt

# GitHub link
github_link = "https://github.com/MichaelTellis/COSMOS_SPARC_project"

# Generate QR code
qr = qrcode.QRCode(
    version=1,
    error_correction=qrcode.constants.ERROR_CORRECT_L,
    box_size=10,
    border=4,
)
qr.add_data(github_link)
qr.make(fit=True)

# Create an image from the QR code
img = qr.make_image(fill_color="black", back_color="white")

# Display the QR code using matplotlib
plt.figure(figsize=(6, 6))
plt.imshow(img, cmap='gray')
plt.axis('off')  # Hide the axes
plt.show()
