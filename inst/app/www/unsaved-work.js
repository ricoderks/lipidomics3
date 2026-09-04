// Warn before the page with the application is left.
//
// Every result of the workflow lives in the memory of the R session that
// belongs to this browser tab, so closing the tab, closing the browser or
// reloading the page throws all of it away. The browser writes the text of the
// warning itself and does not let a page change it.

var lipidomics3UnsavedWork = false;

$(document).ready(function () {
  Shiny.addCustomMessageHandler("lipidomics3-unsaved-work", function (message) {
    lipidomics3UnsavedWork = Boolean(message.unsaved);
  });
});

window.addEventListener("beforeunload", function (event) {
  if (!lipidomics3UnsavedWork) {
    return undefined;
  }

  // Browsers do not agree on how a page asks for the dialog, so all three ways
  // are used.
  event.preventDefault();
  event.returnValue = "";

  return "";
});
