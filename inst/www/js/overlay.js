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

  // Remove Bootstrap modal backdrop so Bayesian progress modal appears without a dark overlay.
  // Use 'shown' (post-animation) so the backdrop exists in DOM when we remove it.
  $(document).on('shown.bs.modal', '.modal', function() {
    $('.modal-backdrop').remove();
  });
  // Bootstrap removes modal-open from body as part of backdrop fadeOut; since we removed
  // the backdrop ourselves, force cleanup so scrolling is restored after modal closes.
  $(document).on('hidden.bs.modal', '.modal', function() {
    $('body').removeClass('modal-open');
  });
});

