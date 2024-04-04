// Function to fetch file structure and populate report.html
function populateIndex() {
  fetchFileStructure('.', function(files) {
    var fileList = document.getElementById('fileList');
    files.forEach(function(file) {
      var listItem = document.createElement('li');
      var link = document.createElement('a');
      link.textContent = file.name;
      link.href = file.path;
      listItem.appendChild(link);
      fileList.appendChild(listItem);
    });
  });
}

// Function to fetch file structure recursively
function fetchFileStructure(path, callback) {
  fetch(path)
    .then(response => response.text())
    .then(text => {
      var parser = new DOMParser();
      var htmlDoc = parser.parseFromString(text, 'text/html');

      // Get immediate folders
      var folders = htmlDoc.querySelectorAll('a[href$="/"]');
      var folderItems = Array.from(folders).map(folder => ({
        name: folder.textContent,
        path: folder.href
      }));

      // Get HTML files
      var links = htmlDoc.querySelectorAll('a[href$=".html"]');
      var fileItems = Array.from(links).map(link => ({
        name: link.textContent,
        path: link.href
      }));

      // Merge and pass the result to the callback
      var files = folderItems.concat(fileItems);
      callback(files);
    })
    .catch(error => {
      console.error('Error fetching file structure:', error);
      callback([]);
    });
}

// Call the function to populate the report.html
populateIndex();
