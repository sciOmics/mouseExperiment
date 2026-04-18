// Mobile device incompatibility overlay
$(document).ready(function() {
  var overlay = $(
    '<div id="mobile-overlay">' +
      '<div class="mobile-icon">\uD83D\uDCBB</div>' +
      '<h2>Desktop Recommended</h2>' +
      '<p>This dashboard requires a larger screen for the best experience. ' +
      'Please visit on a desktop or laptop computer.</p>' +
      '<button class="mobile-continue-btn">Continue anyway</button>' +
    '</div>'
  );
  $('body').append(overlay);

  overlay.on('click', '.mobile-continue-btn', function() {
    overlay.addClass('dismissed');
  });
});
