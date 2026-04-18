// Mobile device incompatibility overlay
$(document).ready(function() {
  var dismissed = false;
  var overlay = $(
    '<div id="mobile-overlay" style="display:none">' +
      '<div class="mobile-icon">\uD83D\uDCBB</div>' +
      '<h2>Desktop Recommended</h2>' +
      '<p>This dashboard requires a larger screen for the best experience. ' +
      'Please visit on a desktop or laptop computer.</p>' +
      '<button class="mobile-continue-btn">Continue anyway</button>' +
    '</div>'
  );
  $('body').append(overlay);

  overlay.on('click', '.mobile-continue-btn', function() {
    dismissed = true;
    overlay.hide();
  });

  var mql = window.matchMedia('(max-width: 768px)');
  function check(e) {
    if (e.matches && !dismissed) {
      overlay.css('display', 'flex');
    } else {
      overlay.hide();
    }
  }
  mql.addEventListener('change', check);
  check(mql);
});
