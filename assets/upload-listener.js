
// Clear expired upload folder from localStorage (after 30 minutes)
const EXPIRY_MINUTES = 2;
const savedUpload = JSON.parse(localStorage.getItem("upload-folder") || "null");

if (savedUpload?.timestamp) {
  const now = Date.now();
  const ageMinutes = (now - savedUpload.timestamp) / (1000 * 60);
  if (ageMinutes > EXPIRY_MINUTES) {
    console.log("🧹 Clearing expired upload-folder from localStorage.");
    localStorage.removeItem("upload-folder");
  }
}

window.addEventListener("message", function (event) {
  if (event.data?.type === "uploadComplete") {
    const folder = event.data.folder;
    console.log("✅ Upload complete. Folder received:", folder);

    const alreadyStored = JSON.parse(localStorage.getItem("upload-folder") || "null");

    if (alreadyStored !== folder) {
      localStorage.setItem("upload-folder", JSON.stringify(folder));
      console.log("💾 Stored new folder in localStorage:", folder);

      // Optional: notify Dash immediately
      const store = document.querySelector('[id="upload-folder"]');
      if (store && store.setProps) {
        store.setProps({ data: folder });
        console.log("📡 Sent folder to Dash via setProps");
      }

      // Force reload to let Dash load from localStorage
      setTimeout(() => {
        console.log("🔄 Reloading page to reflect new upload");
        window.location.reload();
      }, 500);
    } else {
      console.warn("⚠️ Same folder already stored. Skipping reload.");
    }
  }
});
