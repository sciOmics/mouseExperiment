// Mouse background overlay and diagnostics tab disabling
// Extracted from inline JS in app.R (issue 3.3)

$(document).ready(function() {
  var layer = $('<div id="mouse-bg-layer">' +
    '<img id="mouse-bg-img1" class="mouse-bg-img" src="images/mouse_outline1.png" alt="">' +
    '<img id="mouse-bg-img2" class="mouse-bg-img" src="images/mouse_outline2.png" alt="">' +
    '</div>');
  $('body').append(layer);

  // Disable Diagnostics tab — pending re-implementation
  function disableDiagnosticsTab() {
    var $diagLi = $('#tg_results_tabs li').filter(function() {
      return $(this).find('a[data-value="Diagnostics"]').length > 0;
    });
    if ($diagLi.length) {
      $diagLi.addClass('disabled').attr('title', 'Diagnostics — coming in a future release');
      $diagLi.find('a').css({'pointer-events': 'none', 'cursor': 'not-allowed'});
    }
  }
  $(document).one('shiny:sessioninitialized', disableDiagnosticsTab);

  // MutationObserver fires synchronously the moment Bootstrap appends .modal-backdrop
  // to the body — before any CSS transition has a chance to render — and applies
  // inline !important styles (highest possible CSS priority).
  // The backdrop stays in the DOM so Bootstrap's modal-open cleanup runs normally.
  var backdropObserver = new MutationObserver(function(mutations) {
    for (var i = 0; i < mutations.length; i++) {
      var nodes = mutations[i].addedNodes;
      for (var j = 0; j < nodes.length; j++) {
        var node = nodes[j];
        if (node.nodeType === 1 && node.classList && node.classList.contains('modal-backdrop')) {
          node.style.setProperty('background-color', 'transparent', 'important');
          node.style.setProperty('pointer-events', 'none', 'important');
        }
      }
    }
  });
  backdropObserver.observe(document.body, { childList: true });

});

